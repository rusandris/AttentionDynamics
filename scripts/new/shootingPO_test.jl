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


N = 3 
g = 5
w = 0.5

#N = 10 
#g = 100.0
#w = 0.25
wc1 = (N-1)/g

Ttr = 5000
T = 100
#generate standard init cond
p = [g,w,N]
u0 = rand(2N)
fp_guess = ones(2N)
Δt = 1e-2


dsmax = 0.1 #max step size along bif curve
ds = 0.01 #initial step along bif curve
detect_bifurcation = 2 #detection level
max_steps = 50 #max nr of continuation steps

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
x2_ts = [sol.u[i][2] for i in 1:length(sol)]
pl = plot(sol.t[1:2000],x1_ts[1:2000],label="N=$N")
plot!(sol.t[1:2000],x2_ts[1:2000],label="N=$N")


T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷3,ϵ=0.2)
@show T_period

probsh, cish = generate_ci_problem( ShootingProblem(M=200),
	prob, prob_de, sol, T_period; alg = Rodas5(), reltol = 1e-8, abstol = 1e-10)

opts_po_cont = ContinuationPar(p_min = 0.35, p_max = 0.4, ds = -ds,
        newton_options = NewtonPar(;eigsolver=DefaultEig(;which=abs)), save_to_file=false,
        save_eigenvectors=false,save_sol_every_step=100,save_eig_every_step = 1,nev=2N,
        dsmax = dsmax,detect_bifurcation=detect_bifurcation, max_steps = max_steps, tol_stability = 1e-10)

brpo = continuation(probsh, cish, PALC(), opts_po_cont,; 
	verbosity = 3 ,plot = false,
    finalise_solution = fin_sol,funcs_po...)