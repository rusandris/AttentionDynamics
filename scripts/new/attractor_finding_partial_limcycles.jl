using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using Attractors

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.3
p = [g,w,N]

u0 = rand(2N) .* 50
ds = create_ds(u0,p;diffeq=diffeq)

#test trajectory for efficient phase space grid choice
traj,_ = trajectory(ds,500;Ttr=1000,Δt=1e-3)
Plots.plot(traj[:,1],traj[:,2])
#need max finding along each dim
max_values = maximum(Matrix(traj);dims=1) #ugly max finding
grid_size=200
cons_recurrences = 500 #nr of recurrences 
locate_steps = 10000 #nr of steps after detection
max_ensemble = 50 #nr of initial conditions
grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid;Ttr=1000,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps)

mapper(rand(2N) .* 50)