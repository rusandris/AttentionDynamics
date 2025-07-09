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

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-2,adaptive=true,maxiters=1e15)
#N = 10
#g = 18.0
#w = 0.66

N = 5
g = 20.0
w = 0.28
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

traj,t = trajectory(ds,4000,initial_state(ds);Ttr=0,Δt=0.1)
crs = [contraction_rate(traj[i],current_parameters(ds)) for i in 1:length(traj)]


plot(t[1:4000],traj[1:4000,1])

plot(t,traj[:,1])

plot!(crs)

simple_transient_time(ds,T,initial_state(ds); threshold=0.5, ϵ=1e-2, return_ts = false)