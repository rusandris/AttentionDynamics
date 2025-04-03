using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

Ns = [3:6;]
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

w = 0.5
grid_size = 1000 #phase space grid cell size
grid_edges = [0.0,25.0]

T_test = 1000
Ttr=5000
cons_recurrences = 5000; 
locate_steps = 10000;
max_ensemble = 1000
attractor_finding_results = []
#=
grid = (
    range(0.0, 5.0; length = 3), # x1
    range(0.0, 5.0; length = 3), # x2
    range(0.0, 5.0; length = 3), # x3
    range(0.0, 3.0; length = 3), # y1
    range(0.0, 3.0; length = 3), # y2
    range(0.0, 3.0; length = 3), # y3
)
=#

for (n,N) in enumerate(Ns)
	@show N
    g = g_N * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N
    w = wc1 + (wc2-wc1)/2 
    
    @show wc1,wc2,g,w
    
	ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)

    
    #test trajectory for efficient phase space grid choice
    traj,_ = trajectory(ds,T_test;Ttr=Ttr)
    #need max finding along each dim
    max_values = maximum(Matrix(traj);dims=1) #ugly max finding
    
    grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])
    
    
    mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, attractor_locate_steps = locate_steps)
    nr_attractors = 0
    i=0
    ic_ens_idxs = Int[]
    attr_IDs = Int[]
    while (nr_attractors < factorial(N-1) ) || (i > max_ensemble)
    	attr_ID = mapper(rand(2N))
    	nr_attractors = length(mapper.bsn_nfo.attractors)
    	i+=1
    	@show attr_ID
    	push!(ic_ens_idxs,i)
    	push!(attr_IDs,attr_ID)
    end
    push!(attractor_finding_results,(ic_ens_idxs,attr_IDs))
end


