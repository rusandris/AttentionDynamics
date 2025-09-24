using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using Attractors
using Plots
using DelimitedFiles
cd(@__DIR__)
#inspectdr()

savedir = "../../data/attractor_mapper_N5/"
mkpath(savedir)

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 5
g = 6.0

#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N

#original parametrization
K = 1.0;α = 5*1e-2;rc=1.0
competition_strength(w,rc,K,α,N) = w/(1-w) * rc/(α*K*(N-1))

production_rate(g,α) = g*α
rp = production_rate(g, α)
rc = rp

c_crit1 = competition_strength(wc1, rc, K, α, N)
c_crit2 = competition_strength(wc2, rc, K, α, N)
p = [rp,rc,K,4.0,α]

u0 = rand(2N) 
ds = create_original_ds(u0,p;diffeq=diffeq)

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


nr_ens = 10 #nr of initial conditions
attractor_ids = zeros(nr_ens)
init_conds = []
for i in 1:nr_ens
    ic = rand(2N)
    push!(init_conds,ic)
    attr_id = mapper(ic)
    @show i, attr_id
    attractor_ids[i] = attr_id
end

attrs = extract_attractors(mapper)
found_attrs = keys(attrs)
init_conds_idxs_to_attr = [findfirst(i -> i == a,attractor_ids) for a in found_attrs]
u0s = init_conds[init_conds_idxs_to_attr]

writedlm(savedir*"initial_conditions_different_attr_N$N"*"_g$g" * ".txt",stack(u0s))
