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

#---------------------parameters--------------------
Ttr = 5000 #transient time 
T = 400 #saved time

#---------------systems parameters-------------------
N = 5
g = 20.0

#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N


#---------------------------data saving--------------------------

cont_dir = "../../data/bifkit_data_orig_N5_trap_branch2/"
mkpath(cont_dir)

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

#----------define some local parameters------------

#dir for cont result
#cont_dir = "../../bifkit_data_gN5_210/"  
#mkpath(cont_dir)

Δt = 1e-3 #lim cycle sol sampling
dsmax = 0.01 #max step size along bif curve
ds = 0.001 #initial step along bif curve
detect_bifurcation = 2 #detection level
max_steps = 1000 #max nr of continuation steps
M = 650 #discretization for trapezoid method

#nev (nr of eigs) needs to be set every time? to 2N

#--------------------------------continuation for each N-----------------------------
#loop through N values

wc1 = (N-1)/g
wc2 = (N-1)/N


#---------------------------------system parameters------------------------------------

K = 1.0;α = 5*1e-2;rc=1.0
competition_strength(w,rc,K,α,N) = w/(1-w) * rc/(α*K*(N-1))
c_crit1 = competition_strength(wc1, rc, K, α, N)
c_crit2 = competition_strength(wc2, rc, K, α, N)

global_params = [N, g, K, α, c_crit1,c_crit2]
#save params that apply to every system size
writedlm(cont_dir * "global_params_Ns_g.txt", global_params)


production_rate(g,α) = g*α
rp = production_rate(20.0, α)
rc = rp
p_orig = [rp,rc,K,1.1,α]

c_min = c_crit1 
c_max = c_crit2

f1,f2 = rescaling_factors(rp, rc, K, 1.1, α, N)


#-----------------coexistence fixed point continuation-----------------
#follow coex fp up until the loss of stability (high d Hopf)
fp_guess = f1 .* ones(2N) #coexistence fp     

#def problem 
prob = BifurcationProblem(AD_orig,
        fp_guess, # initial condition x0
        p_orig, # set of parameters
        4,     # parameter index for continuation
        record_from_solution = (x, p; k...) -> x[1])

#p_min = 0.1, p_max = wc1+0.01
opts_fp_cont = ContinuationPar(p_min = 1.0, p_max = 1.3, ds = ds,
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

#------------------solve for a lim cycle---------------
#take lim cycle closer to wc2


c_cent = 6.0 

#params for continuation
p_orig = p_orig = [rp,rc,K,c_cent,α]
u0 = rand(2N)

#transient
prob_de = ODEProblem(att_din_original!, u0, (0,Ttr), p_orig;saveat=[Ttr])
sol = DifferentialEquations.solve(prob_de, Rodas5())

#higher precision, shorter traj (couple of periods)
prob_de = ODEProblem(att_din_original!, sol.u[end], (0,T), p_orig, reltol = 1e-8, abstol = 1e-10)
sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)
#display(plot(sol))

#approximate period 
x1_ts = [sol.u[i][1] for i in 1:length(sol)]
display(plot(sol.t,x1_ts,label="N=$N"))
T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷2,ϵ=0.5)
@show T_period

#setup for periodic orbit (lim cycle) continuation

#trapezoid method continuation 
prob_orbit, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M=M), prob, sol, T_period)

#collocation method continuation 
#prob_orbit, orbit_guess = generate_ci_problem(PeriodicOrbitOCollProblem(30, 5; meshadapt=true),
#    prob, sol, T_period)

#save max of x1
#perform continuation w_cent -> w_min direction (backwards)
#------periodic orbit continuation parameters---------
opts_po_cont = ContinuationPar(p_min=c_min, p_max=c_cent, ds=-ds,
    newton_options=NewtonPar(; eigsolver=DefaultEig(; which=abs)), save_to_file=false,
    save_eigenvectors=false, save_sol_every_step=100, save_eig_every_step=1, nev=2N,
    dsmax=dsmax, detect_bifurcation=detect_bifurcation, max_steps=max_steps, tol_stability=1e-10)

brpo1 = continuation(prob_orbit, orbit_guess, PALC(), opts_po_cont, ;
verbosity = 3,filename=cont_dir*"limit_cycle_branch_N$N" ,plot = false,funcs_po...)

#perform continuation w_cent -> w_max direction (forwards)
opts_po_cont = ContinuationPar(p_min=c_cent, p_max=c_max, ds=ds,
    newton_options=NewtonPar(; eigsolver=DefaultEig(; which=abs)), save_to_file=false,
    save_eigenvectors=false, save_sol_every_step=100, save_eig_every_step=1, nev=2N,
    dsmax=dsmax, detect_bifurcation=detect_bifurcation, max_steps=max_steps, tol_stability=1e-10)

brpo2 = continuation(prob_orbit, orbit_guess, PALC(), opts_po_cont, ;
    verbosity=3, filename=cont_dir * "limit_cycle_branch_N$N", plot=false, funcs_po...)

#-------------------extract data from continuation res-------------------------

for (i,brpo) in enumerate([brpo2])
    #solution info: maxs,mins of x1 of the solution for each param 
    params = brpo.branch.param 
    mins = brpo.branch.min
    maxs = brpo.branch.max
    periods = brpo.branch.period
    stability = brpo.stable 

    writedlm(cont_dir*"contres_po_sols$(i)_N$N"*".txt",hcat(params,mins,maxs,periods,stability))



    #bifurcation point info:
    bifpoints = []
    for bp in brpo.specialpoint
        push!(bifpoints,[bp.type,bp.param,bp.idx])
    end
    writedlm(cont_dir*"contres_po_bifp$(i)_N$N" * ".txt", bifpoints)
end


#eig info: max of real parts of floquet exps
#floq_exps = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
#floq_mult = [exp.(floq_exps[i] .* periods[i]) for i in 1:length(floq_exps)]

#=
#filter out zero exponents
floq_exps_no_zero = [filter(μ -> !(iszero_approx(μ)),floq_exps[i]) for i in 1:length(floq_exps)]

@show length(floq_exps_no_zero)

real_floq_exps = []
imag_floq_exps = []
for i in 1:length(floq_exps_no_zero)
    if length(floq_exps_no_zero[i]) != 0
        push!(real_floq_exps,real(floq_exps_no_zero[i]))
        push!(imag_floq_exps,imag(floq_exps_no_zero[i]))
    end
end

maxreal_floq_exps = maximum.(real_floq_exps)
maximag_floq_exps = minimum.(imag_floq_exps)
=#

#@show length.(floq_exps)
#writedlm(cont_dir*"contres_po_eig_maxfloq_N$N"*".txt",hcat(maxreal_floq_exps,maximag_floq_exps))
#writedlm(cont_dir*"contres_po_eig_allfloq_N$N"*".txt",hcat(stack(real.(floq_exps))',stack(imag.(floq_exps))'))
