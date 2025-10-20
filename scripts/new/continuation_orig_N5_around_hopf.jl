using Revise
using BifurcationKit
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using Plots
#using DynamicalSystems
using ChaosTools
cd(@__DIR__)

#--------TODO-------------
#how to plot bifurcation diagrams of multiple limit cycles? 
#save bifurcation diagrams of many distinct limit cycles 
#save an initial condition mapper (3-4 attractors)
#save continuation results of each initial condition
#plot together 

#---------------------parameters--------------------
Ttr = 5000 #transient time 
T = 2500 #saved time

#---------------systems parameters-------------------
N = 5
g = 6.0

#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N


#----------------------------helper functions-----------------------------

isunity(z) = isapprox(real(z),1.0;atol=1e-1) && isapprox(imag(z),0.0;atol=1e-1)
iszero_approx(z) = isapprox(real(z),0.0;atol=1e-3) && isapprox(imag(z),0.0;atol=1e-3)

#---------------------define system-----------------------
#general so dimension can be changed

function AD_orig(u,p)
    dim = length(u)
    rp,rc,K,c,α = p
    N = Int(length(u)/2)
    du = zeros(Real,dim)
    for i in 1:2:dim
        du[i] = rp*u[i]*(1 - rc/K*u[i+1] - c*skipsum(u,i))
        du[i+1] = u[i] - α*u[i+1]
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

function lim_cycle_continuation(prob,sol,p_min,p_max; M = 200, ds = 0.01, max_steps=500,ϵ = 0.4,detect_bifurcation=2)

    dim = length(sol.u[1])
    #approximate period 
    x1_ts = [sol.u[i][1] for i in 1:length(sol)]
    display(plot(sol.t,x1_ts,label="N=$N"))
    T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷2,ϵ=ϵ)
    @show T_period

    #setup for periodic orbit (lim cycle) continuation

    #trapezoid method continuation 
    prob_orbit, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M=M), prob, sol, T_period)

    #collocation method continuation 
    #prob_orbit, orbit_guess = generate_ci_problem(PeriodicOrbitOCollProblem(30, 5; meshadapt=true),
    #    prob, sol, T_period)

    #save max of x1
    #perform continuation max -> min direction (backwards)
    #------periodic orbit continuation parameters---------
    opts_po_cont = ContinuationPar(p_min=c_min, p_max=c_max, ds=-ds,
        newton_options=NewtonPar(; eigsolver=DefaultEig(; which=abs)), save_to_file=false,
        save_eigenvectors=false, save_sol_every_step=100, save_eig_every_step=1, nev=dim,
        dsmax=dsmax, detect_bifurcation=detect_bifurcation, max_steps=max_steps, tol_stability=1e-10)

    brpo = continuation(prob_orbit, orbit_guess, PALC(), opts_po_cont, ;
    verbosity = 3,plot = false,funcs_po...)

    return brpo


end


function extract_continuation_result(brpo;cont_dir)
    #-------------------extract data from continuation res-------------------------

    #solution info: maxs,mins of x1 of the solution for each param 
    params = brpo.branch.param 
    mins = brpo.branch.min
    maxs = brpo.branch.max
    periods = brpo.branch.period
    stability = brpo.stable 

    writedlm(cont_dir*"contres_po_sols_N$N"*".txt",hcat(params,mins,maxs,periods,stability))

    #bifurcation point info:
    bifpoints = []
    for bp in brpo.specialpoint
        push!(bifpoints,[bp.type,bp.param,bp.idx])
    end
    writedlm(cont_dir*"contres_po_bifp_N$N" * ".txt", bifpoints)

end

#----------define some local parameters------------

#dir for cont result
#cont_dir = "../../bifkit_data_gN5_210/"  
#mkpath(cont_dir)

Δt = 1e-2 #lim cycle sol sampling
dsmax = 0.01 #max step size along bif curve
ds = 0.001 #initial step along bif curve
detect_bifurcation = 2 #detection level
max_steps = 300 #max nr of continuation steps
M = 450 #discretization for trapezoid method

#nev (nr of eigs) needs to be set every time? to 2N

#--------------------------------continuation for each N-----------------------------
#loop through N values

wc1 = (N-1)/g
wc2 = (N-1)/N

#---------------------------------system parameters------------------------------------

K = 1.0;α = 5*1e-2;rc=1.0
competition_strength(w,rc,K,α,N) = w/(1-w) * rc/(α*K*(N-1))

production_rate(g,α) = g*α
rp = production_rate(g, α)
rc = rp

c_crit1 = competition_strength(wc1, rc, K, α, N)
c_crit2 = competition_strength(wc2, rc, K, α, N)

#------------------solve for a lim cycle---------------
    #take lim cycle closer to wc2

c_min = c_crit1 - 0.3
c_max = c_crit1 + 1.0
#c_max = c_crit1 + 4.0 
p_orig = [rp,rc,K,c_max,α]

f1,f2 = rescaling_factors(rp, rc, K, c_max, α, N)

fp_guess = f1 .* ones(2N) #coexistence fp     

#def problem 
prob = BifurcationProblem(AD_orig,
        fp_guess, # initial condition x0
        p_orig, # set of parameters
        4,     # parameter index for continuation
        record_from_solution = (x, p; k...) -> x[1])

#read initial conditions 
dir_ic = "../../data/attractor_mapper_N5/" 
u0s_matrix = readdlm(dir_ic*"initial_conditions_different_attr_N5_g6.0.txt")
u0s = [u0s_matrix[:,i] for i in 1:size(u0s_matrix)[2] ]

#data saving
dir_cont = "../../data/bifkit_data_orig_N5_around_hopf/"
mkpath(dir_cont)

global_params = [N, g, K, α, c_crit1,c_crit2]
#save params that apply to every system size

#loop through initial conditions
for (i,u0) in enumerate(u0s)  

    @show i

    #transient
    prob_de = ODEProblem(att_din_original!, u0, (0,Ttr), p_orig;saveat=[Ttr])
    sol = DifferentialEquations.solve(prob_de, Rodas5())

    #higher precision, shorter traj (couple of periods)
    prob_de = ODEProblem(att_din_original!, sol.u[end], (0,T), p_orig, reltol = 1e-8, abstol = 1e-10)
    sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)
    #display(plot(sol))

    brpo = lim_cycle_continuation(prob,sol,c_min,c_max; M = M, ds = ds, max_steps=max_steps,ϵ = 0.4,detect_bifurcation=2)
    dir_cont_i = dir_cont*"attr$i"*"/"
    mkpath(dir_cont_i)
    extract_continuation_result(brpo;cont_dir=dir_cont_i)
    writedlm(dir_cont_i * "global_params_Ns_g.txt", global_params)
end



#=
#-----------------coexistence fixed point continuation-----------------
#follow coex fp up until the loss of stability (high d Hopf)


c_min_fp = c_crit1-1.0
c_max_fp = c_crit2 + 1.0

p_orig = [rp,rc,K,c_min_fp,α]

f1,f2 = rescaling_factors(rp, rc, K, c_min_fp, α, N)

fp_guess = f1 .* ones(2N) #coexistence fp     

#def problem 
prob = BifurcationProblem(AD_orig,
        fp_guess, # initial condition x0
        p_orig, # set of parameters
        4,     # parameter index for continuation
        record_from_solution = (x, p; k...) -> x[1])

#p_min = 0.1, p_max = wc1+0.01
opts_fp_cont = ContinuationPar(p_min = c_min_fp, p_max = c_max_fp, ds = ds,
    newton_options = NewtonPar(;eigsolver=DefaultEig(;which=real)), save_to_file=false,
    save_eigenvectors=false,save_sol_every_step=100,save_eig_every_step = 1,
    dsmax = dsmax,detect_bifurcation=detect_bifurcation, max_steps = max_steps, tol_stability = 1e-10)

brfp = continuation(prob, PALC(), opts_fp_cont; 
verbosity = 3 ,plot = false)

#-------------------extract data from continuation res-------------------------
#solution info: x1 of the solution for each param 
params = brfp.branch.param 
x = brfp.branch.x
stability = brfp.stable 

writedlm(cont_dir*"contres_fp_sols_N$N"*".txt",hcat(params,x,stability))

#bifurcation point info:
bifpoints = []
for bp in brfp.specialpoint
    push!(bifpoints,[bp.type,bp.param,bp.idx])
end
writedlm(cont_dir*"contres_fp_bifp_N$N"*".txt",bifpoints)

=#