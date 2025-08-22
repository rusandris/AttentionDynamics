using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using Attractors
using Plots
using DelimitedFiles
#inspectdr()

savedir = "../../data/"
mkpath(savedir)

#diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.3
p = [g,w,N]

u0 = rand(2N) .* 50
ds = create_ds(u0,p;diffeq=diffeq)

#test trajectory for efficient phase space grid choice
traj,t = trajectory(ds,500;Ttr=1000,Δt=1e-2)

#need max finding along each dim
min,max = extrema(Matrix(traj)) #ugly max finding

#--------------grid parameters------------------
grid_size = 200
cons_recurrences = 500 #nr of recurrences 
locate_steps = 10000 #nr of steps after detection
Ttr = 2000

grid = Tuple([range(min,max;length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps)


nr_ens = 1000 #nr of initial conditions
attractor_ids = zeros(nr_ens)
for i in 1:nr_ens
  attr_id = mapper(rand(2N) .* 50)
  @show i, attr_id
  attractor_ids[i] = attr_id
end

attrs = extract_attractors(mapper)
nr_actives = zeros(nr_ens)

for i in keys(attrs)
    #count nr of active variables in limit cycle
    A = Matrix(attrs[i])
    nr_active = count(maximum(A[:,1:2:end];dims=1) .> 1.0)
    nr_actives[i] = nr_active
end

writedlm(savedir*"attractor_finding_results_partial_limcycles_Rodas5_solver_N$N"*"_g$g" * "cons_rec500" * ".txt",hcat(1:nr_ens,attractor_ids,nr_actives))
