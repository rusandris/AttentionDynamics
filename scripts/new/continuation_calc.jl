using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
#using DynamicalSystems
using ChaosTools
cd(@__DIR__)

#include standard parameters and init conditions
#run necessary mkpaths
include("bif_floq_od_method_parameters.jl")


isunity(z) = isapprox(real(z),1.0;atol=1e-1) && isapprox(imag(z),0.0;atol=1e-1)
iszero_approx(z) = isapprox(real(z),0.0;atol=1e-2) && isapprox(imag(z),0.0;atol=1e-2)

#---------------------define system-----------------------
#general so dimension can be changed

function skipsum(a,i)
	x = a[1:2:end]
    return sum(x) - a[i]
end

function F_ad(u,p)
    dim = length(u)
    N = dim/2
    du = Real[]
    for i in 1:2:(dim)
        dxi = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(N-1)*skipsum(u,i))
        dyi = u[i] - u[i+1]
        push!(du,dxi,dyi)
    end
    du
end

#-------------------helper recording and plotter functions---------------------

funcs_po = (record_from_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, p.p))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
		# plot!(br; subplot = 1, putspecialptlegend = false)
	end)

#----------define some local parameters------------

Δt = 1e-2 #lim cycle sol sampling
dsmax = 0.02 #max step size along bif curve
ds = -0.02 #initial step along bif curve
detect_bifurcation = 2 #detection level
max_steps = 250 #max nr of continuation steps
M = 250 #discretization for trapezoid method

#nev (nr of eigs) needs to be set every time? to 2N

#--------------------------------continuation for each N-----------------------------
#loop through N values

for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #generate standard init cond
    u0 = generate_initial_condition(N)

    #params for continuation
    w_min = w_min_ratio * wc1
    w_max = w_max_ratio * wc2
    p = [g,w_max,N]

    #def problem 
    prob = BifurcationProblem(F_ad,
            u0, # initial condition x0
            p, # set of parameters
            2;     # parameter index for continuation
            )

    #------------------solve for a lim cycle---------------
    #take lim cycle closer to wc2

    #transient
    prob_de = ODEProblem(att_din!, u0, (0,Ttr), p;saveat=[Ttr])
    sol = DifferentialEquations.solve(prob_de, Rodas5())

    #higher precision, shorter traj (couple of periods)
    prob_de = ODEProblem(att_din!, sol.u[end], (0,100.), p, reltol = 1e-8, abstol = 1e-10)
    sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)
    #display(plot(sol))

    #approximate period 
    x1_ts = [sol.u[i][1] for i in 1:length(sol)]
    display(plot(sol.t,x1_ts,label="N=$N"))
    T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷3,ϵ=0.2)
    @show T_period

    #setup for periodic orbit (lim cycle) continuation
    #------periodic orbit continuation parameters---------
    opts_po_cont = ContinuationPar(p_min = w_min, p_max = w_max, ds = ds,
        newton_options = NewtonPar(;eigsolver=DefaultEig(;which=abs)), save_to_file=false,
        save_eigenvectors=false,save_sol_every_step=100,save_eig_every_step = 1,
        dsmax = dsmax,detect_bifurcation=detect_bifurcation, max_steps = max_steps, tol_stability = 1e-10)
    
    #trapezoid method continuation 
    probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = M),prob, sol, T_period)
    
    #save max of x1
    #perform continuation wc_max -> w_min direction (backwards)
    brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont,; 
	verbosity = 3,filename=cont_dir*"limit_cycle_branch_N$N" ,plot = false,funcs_po...)

    #-------------------extract data from continuation res-------------------------

    #solution info: maxs,mins of x1 of the solution for each param 
    params = brpo.branch.param 
    mins = brpo.branch.min
    maxs = brpo.branch.max
    periods = brpo.branch.period
    stability = brpo.stable 

    writedlm(cont_dir*"contres_sols_N$N"*".txt",hcat(params,mins,maxs,periods,stability))

    #eig info: max of real parts of floquet exps
    floq_exps = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
    floq_mult = [exp.(floq_exps[i].*periods[i]) for i in 1:length(floq_exps)]

    #filter out zero exponents
    floq_exps_no_zero = [filter(μ -> !(iszero_approx(μ)),floq_exps[i]) for i in 1:length(floq_exps)]

    maxreal_floq_exps = maximum.(real.(floq_exps_no_zero))
    writedlm(cont_dir*"contres_eig_maxreal_N$N"*".txt",maxreal_floq_exps)

    #bifurcation point info:
    bifpoints = []
    for bp in brpo.specialpoint
        push!(bifpoints,[bp.type,bp.param,bp.idx])
    end
    writedlm(cont_dir*"contres_bifp_N$N"*".txt",bifpoints)


    #=
    #should filter out unity values
    @show length(floq_mult)
    floq_mult_no_unity = [filter(λ -> !(isunity(λ)),floq_mult[i]) for i in 1:length(floq_mult)]

    real_parts = stack([real.(floq_mult[i]) for i in 1:length(floq_mult)])'
    imag_parts = stack([imag.(floq_mult[i]) for i in 1:length(floq_mult)])'

    push!(imag_parts_allN,imag_parts)
    push!(real_parts_allN,real_parts)
    =#
   
end


