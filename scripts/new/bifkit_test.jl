using Revise
using BifurcationKit
using Plots
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra

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
=#

#use DiffEq to solve for a periodic trajectory

N = 5 
g = 7.5
w = 0.75
wc1 = (N-1)/g
wc2 = (N-1)/N
p = [g,w,N]
u0 = rand(2N)
Δt = 1e-2

prob = BifurcationProblem(F_ad,
        u0, # initial condition x0
        p, # set of parameters
        2;     # parameter index for continuation
        record_from_solution = (x,p; k...) -> x[1])

#br = continuation(prob, PALC(), ContinuationPar(p_min = 0.5, p_max = 0.8))


#transient
prob_de = ODEProblem(att_din!, u0, (0,1000.0), p;saveat=[1000.0])
sol = DifferentialEquations.solve(prob_de, Rodas5())

#higher precision, shorter traj (couple of periods)
prob_de = ODEProblem(att_din!, sol.u[end], (0,20.), p, reltol = 1e-8, abstol = 1e-10)
sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)

plot(sol)
#approximate period 
T_period = 13.0

probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = 300),prob, sol, T_period)
opts_po_cont = ContinuationPar(p_min = 0.55, p_max = 0.75, ds = -0.02, newton_options=NewtonPar(;eigsolver=DefaultEig(;which=imag)), 
    dsmax = 0.02, detect_bifurcation=1, nev = 3, max_steps = 20, tol_stability = 1e-10)

brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont; 
	verbosity = 3, plot = true,finalise_solution = (z, tau, step, contResult; k...) -> begin
    @info "Floquet exponents:"
    (Base.display(contResult.eig[end].eigenvals) ;true)
    end)

periods = brpo.branch.period 
floqs = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
floq_mult = [exp.(floqs[i].*periods[i]) for i in 1:length(floqs)]
floq_mult_max_abs = [maximum(abs.(exp.(floqs[i].*periods[i]))) for i in 1:length(floqs)]

#plot eigenvalues in the complex plane
#start and end of continuation
real_parts_start = real.(floq_mult[1])
imag_parts_start = imag.(floq_mult[1])

real_parts_end = real.(floq_mult[end])
imag_parts_end = imag.(floq_mult[end])

plot(real_parts_start,imag_parts_start,st=:scatter,mc=:gray10,xlims=(-2,2),ylims=(-2,2))
plot!(real_parts_end,imag_parts_end,st=:scatter,mc=:orange,xlims=(-2,2),ylims=(-2,2))

#extract w params
ws = brpo.branch.param

real_parts = stack([real.(floq_mult[i]) for i in 1:length(floq_mult)])'
imag_parts = stack([imag.(floq_mult[i]) for i in 1:length(floq_mult)])'
plot(ws,imag_parts[:,1:3],st=:scatter,ms=2,markerstrokewidth=0.0)
plot(ws,real_parts[:,1:3],st=:scatter,ms=2,markerstrokewidth=0.0)
abs_value = sqrt.(real_parts[:,2] .^2 + imag_parts[:,2] .^2)
plot(ws,abs_value,st=:scatter,ms=2,markerstrokewidth=0.0)

#=
#test eigenvals of the Jacobian of the cont system
#loop through the whole lim cycle
eigvals_all_jac = []
for i in 1:length(sol) 
    u = sol[i]
    J = build_jacobian(u,p)
    eigvals_jac = eigvals(J)
    push!(eigvals_all_jac,eigvals_jac)
end
=#