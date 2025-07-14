using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")
include("transient_functions.jl")


#old comments for parameter setting
#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients


#how to test if the trajectories reached the attractor?
#previous method of finding regular time intervals between peaks isn't reliable
#1.periodically calc lyap exp? (this might fail if the system has long intermittency)
#-> test: plot LLE along trajectory to see if it fluctuates
#-> second test: plot phase space contraction rate to see if it fluctuates
#-> 2.test recurrence (on poincare plane)

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.66

N = 5
g = 20.0
w = 0.21
u0 = rand(2N)

ds = create_ds(u0,[g,w,N];diffeq=diffeq)
T = 500
Ttr = 1500

plane = (2,1.0)
pmap = PoincareMap(ds,plane) 

xs = []
for _ in 1:6500
    step!(pmap)
    push!(xs,current_state(pmap)[1])
end
plot(xs,st=:scatter,ms=3)

traj,t = trajectory(ds,200,initial_state(ds);Ttr=0,Δt=0.001)
crs = [contraction_rate(traj[i],current_parameters(ds)) for i in 1:length(traj)]
simple_transient_time(ds,T,initial_state(ds); threshold=1.5, period_tol=1e-4, max_tol=Inf, return_ts = false)

plot(t[:],traj[:,1])
plot!(xlims=(600,800))

lyapunov(ds,500;u0=initial_state(ds),Ttr=1000)

plot(t,traj[:,1])

plot!(crs)

lifetimes = []
for i in 1:100
    @show i
    lt = simple_transient_time(ds,T,rand(2N); threshold=1.5, ϵ=1e-4, return_ts = false)
    push!(lifetimes,lt)
end