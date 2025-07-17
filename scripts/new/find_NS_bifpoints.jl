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
iszero_approx(z) = isapprox(real(z),0.0;atol=1e-3) && isapprox(imag(z),0.0;atol=1e-3)

#---------------------define system-----------------------
#general so dimension can be changed

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

#finalizer
function fin_sol(z, tau, step, br; k...)
    if length(br.specialpoint) > 0 
        return br.specialpoint[end].type != :ns
    else
        return true
    end
end
#----------define some local parameters------------

#dir for cont result
cont_dir = "../../data/NS_bifpoints_data_gN10/"  
mkpath(cont_dir)

Δt = 1e-3 #lim cycle sol sampling
dsmax = 0.1 #max step size along bif curve
ds = 0.01 #initial step along bif curve
detect_bifurcation = 2 #detection level
max_steps = 50 #max nr of continuation steps
M = 600 #discretization for trapezoid method
Ttr=1000

#nev (nr of eigs) needs to be set every time? to 2N

N = 10 #system size 
phase_diag_params = readdlm("../../data/phasediag_data_N10_zoomed_test/phasediag_params.txt") #read from phase diagram params
gs = phase_diag_params[end:-1:1,2]
#gs = range(2.0,3.0;length=50) .* (N-1)
gs_norm = gs ./ (N-1) 
#gs = gs[7:end]
#s = range(1.2,2.0;length=40) .* (N-1)

g_critical_curve(w,N) = (N-1)/w
fp_guess = ones(2N) #coexistence fp


#loop over g values
for (i,g) in enumerate(gs)

    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    #g = gval_phase_diag * (N-1)
    #g = g_rel*(N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

      #------------------solve for a lim cycle---------------
    #take lim cycle closer to wc2

    #params for continuation
    w_min = wc1 + 0.1
    #w_max = w_max_ratio * wc2
    #w_min = 0.201
    w_max = 0.45
    @show w_min,w_max

    #generate standard init cond
    p = [g,w_max,N]
    u0 = generate_initial_condition(N)

    #def problem           
    prob = BifurcationProblem(F_ad,
            fp_guess, # initial condition x0
            p, # set of parameters
            2,     # parameter index for continuation
            record_from_solution = (x, p; k...) -> x[1])

    #transient
    prob_de = ODEProblem(att_din!, u0, (0,Ttr), p;saveat=[Ttr])
    sol = DifferentialEquations.solve(prob_de, Rodas5())

    #higher precision, shorter traj (couple of periods)
    prob_de = ODEProblem(att_din!, sol.u[end], (0,200.), p, reltol = 1e-8, abstol = 1e-10)
    sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)
    #display(plot(sol))

    #approximate period 
    x1_ts = [sol.u[i][1] for i in 1:length(sol)]
    #display(plot(sol.t,x1_ts,label="N=$N"))
    T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷3,ϵ=0.2)
    @show T_period

    #setup for periodic orbit (lim cycle) continuation
    #------periodic orbit continuation parameters---------
    opts_po_cont = ContinuationPar(p_min = w_min, p_max = w_max, ds = -ds,
        newton_options = NewtonPar(;eigsolver=DefaultEig(;which=abs)), save_to_file=false,
        save_eigenvectors=false,save_sol_every_step=100,save_eig_every_step = 1,nev=2N,
        dsmax = dsmax,detect_bifurcation=detect_bifurcation, max_steps = max_steps, tol_stability = 1e-10)
    
    #trapezoid method continuation 
    probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = M),prob, sol, T_period)
    
    #save max of x1
    #perform continuation wc_max -> w_min direction (backwards)
    brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont,; 
	verbosity = 3,filename=cont_dir*"limit_cycle_branch_N$N" ,plot = false,
    finalise_solution = fin_sol,funcs_po...)

    #-------------------extract data from continuation res-------------------------

    #solution info: maxs,mins of x1 of the solution for each param 
    params = brpo.branch.param 
    mins = brpo.branch.min
    maxs = brpo.branch.max
    periods = brpo.branch.period
    stability = brpo.stable 

    #bifurcation point info:
    #=
    bifpoints = []
    for bp in brpo.specialpoint
        push!(bifpoints,[bp.type,bp.param,bp.idx])
    end
    =#
    @show brpo.specialpoint[1].type
    open(cont_dir*"contres_po_bifp_ns"*".txt", "a") do io
        writedlm(io, [brpo.specialpoint[1].type brpo.specialpoint[1].param])
    end
  



end