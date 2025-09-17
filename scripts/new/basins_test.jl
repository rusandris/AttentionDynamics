using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using Attractors
using Plots
using DelimitedFiles
cd(@__DIR__)
#inspectdr()

savedir = "../../data/"
mkpath(savedir)

#define parameters, system size

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 3
g = 5.0
w = 0.6
p = [g,w,N]
wc1 = (N-1)/g
wc2 = (N-1)/N

#create ds
u0 = rand(2N) 
ds = create_ds(u0,p;diffeq=diffeq)

#create ds projection
projection = [1,3]
complete_state = zeros(4) #[0.0,0.0,1.0,0.0] #u0[[2,4,5,6]] #ones(2N-2) 
ds_x1x2 = ProjectedDynamicalSystem(ds,projection,complete_state)

#test trajectory for efficient phase space grid choice
traj,t = trajectory(ds,500;Ttr=1000,Δt=1e-2)

#need max finding along each dim
min,max = extrema(Matrix(traj)) #ugly max finding

#--------------grid parameters------------------
grid_size = 20
cons_recurrences = 500 #nr of recurrences 
locate_steps = 1000 #nr of steps after detection
Ttr = 200

grid_attractors = Tuple([range(min,max;length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid_attractors;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps,Δt=0.1)

#define grid for basins 

x1_grid = x2_grid = range(0.01, 2.0; length = 100)
grid_for_basin = (x1_grid,x2_grid)  

basins = basins_of_attraction(mapper, grid_for_basin)