using Revise
using BifurcationKit
using Plots
using AttentionDynamics
using DifferentialEquations

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

#=
N = 5 
g = 7.5
w = 0.5
wc1 = (N-1)/g
wc2 = (N-1)/N

p = [g,w,N]
#u0 = rand(2N)
u0 = ones(2N)

prob = BifurcationProblem(F_ad,
        u0, # initial condition x0
        p, # set of parameters
        2;     # parameter index for continuation
        record_from_solution = (x,p; k...) -> x[1])

br = continuation(prob, PALC(), ContinuationPar(p_min = 0.5, p_max = 0.6))

=#
#use DiffEq to solve for a periodic trajectory

N = 5 
g = 7.5
w = 0.7
wc1 = (N-1)/g
wc2 = (N-1)/N
p = [g,w,N]
u0 = rand(2N)

#transient
prob_de = ODEProblem(att_din!, u0, (0,1000.0), p;saveat=[1000.0])
sol = DifferentialEquations.solve(prob_de, Rodas5())

#higher precision, shorter traj (couple of periods)
prob_de = ODEProblem(att_din!, sol.u[end], (0,20.), p, reltol = 1e-8, abstol = 1e-10)
sol = DifferentialEquations.solve(prob_de, Rodas5())

plot(sol)
#approximate period 
T_period = 10.0

probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = 300),prob, sol, T_period)
opts_po_cont = ContinuationPar(p_min = 0.55, p_max = 0.75, ds = -0.002, dsmax = 0.01, n_inversion = 6, nev = 4, max_steps = 50, tol_stability = 1e-8)

brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont,; 
	verbosity = 3, plot = true,finalise_solution = (z, tau, step, contResult; k...) -> begin
    @info "Floquet exponents:"
    (Base.display(contResult.eig[end].eigenvals) ;true)
    end)

periods = brpo.branch.period 
floqs = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
floq_mult_max_abs = [maximum(abs.(exp.(floqs[i].*periods[i]))) for i in 1:length(floqs)]