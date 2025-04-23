using Revise
using BifurcationKit
using Plots
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using DynamicalSystems
cd(@__DIR__)

isunity(z) = isapprox(real(z),1.0;atol=1e-2) && isapprox(imag(z),0.0;atol=1e-2)

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

#----------define parameters------------
Ns = [3:7;] #nr of competitors dim = 2N 

#select g value from the phase diagram
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
gval_phase_diag = gs[30]
Δt = 1e-2 #lim cycle sol sampling
real_parts_allN = [] #store floq multipliers real parts
imag_parts_allN = [] #store floq multipliers imag parts
special_points = []
ws_step_values = []

#loop through N values
for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N
    #w = wc1 + (wc2-wc1)/2
    #@show wc1,wc2,w,g
    u0 = rand(2N)

    #params for continuation
    w_min = 1.05 * wc1
    w_max = 0.85 * wc2
    p = [g,w_max,N]

    #def problem 
    prob = BifurcationProblem(F_ad,
            u0, # initial condition x0
            p, # set of parameters
            2;     # parameter index for continuation
            record_from_solution = (x,p; k...) -> x[1])

    #------------------solve for a lim cycle---------------
    #take lim cycle closer to wc2

    #transient
    prob_de = ODEProblem(att_din!, u0, (0,5000.0), p;saveat=[5000.0])
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
    #perform continuation wc_max -> w_min direction (backwards)
    probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = 250),prob, sol, T_period)
    opts_po_cont = ContinuationPar(p_min = w_min, p_max = w_max, ds = -0.02,newton_options=NewtonPar(;eigsolver=DefaultEig(;which=abs)),
    dsmax = 0.02,detect_bifurcation=2, nev = 2N, max_steps = 150, tol_stability = 1e-10)

    brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont,; 
	verbosity = 3, plot = false)
    #finalise_solution = (z, tau, step, contResult; k...) -> begin
    #@info "Floquet exponents:"
    #(Base.display(contResult.eig[end].eigenvals) ;true)
    #end)

    push!(special_points,brpo.specialpoint)
    #extract data from continuation branch
    periods = brpo.branch.period 
    @show periods[1]
    floqs = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
    floq_mult = [exp.(floqs[i].*periods[i]) for i in 1:length(floqs)]

    #extract w params for all continuation steps
    ws = brpo.branch.param
    push!(ws_step_values,ws)
    #should filter out unity values
    @show length(floq_mult)
    floq_mult_no_unity = [filter(λ -> !(isunity(λ)),floq_mult[i]) for i in 1:length(floq_mult)]

    real_parts = stack([real.(floq_mult_no_unity[i]) for i in 1:length(floq_mult_no_unity)])'
    imag_parts = stack([imag.(floq_mult_no_unity[i]) for i in 1:length(floq_mult_no_unity)])'

    push!(imag_parts_allN,imag_parts)
    push!(real_parts_allN,real_parts)
end

#------------------------plot Floquet multipliers--------------------------
#plot all non-unity floquet multis (real and imag parts)    
plots_imag = []
plots_real = []
for (i,N) in enumerate(Ns)
    pl_imag = plot(ws_step_values[i],imag_parts_allN[i],mc=:gray10,st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    pl_real = plot(ws_step_values[i],real_parts_allN[i],mc=:gray10,st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    push!(plots_imag,pl_imag)
    push!(plots_real,pl_real)
end
plot(plots_real...,layout=(1,length(Ns)),size=(1200,300),ylims=[-1,1])
plot(plots_imag...,layout=(1,length(Ns)),size=(1200,300),ylims=[-1,1])