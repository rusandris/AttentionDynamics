using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using Plots,LaTeXStrings,GLMakie
plotly()
cd(@__DIR__)
include("../../plotting/plot_functions/orbit_diagrams.jl")
include("../../plotting/defaults.jl")
include("bif_floq_od_method_parameters.jl")

#parameters from phase diagram 
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

#--------------------parameters for attractor finding via recurrences---------------
grid_size = 200 #phase space grid cell size
T_test = 1000 #test trajectory for grid edges 
#Ttr = 10000 #transient time
cons_recurrences = 500 #nr of recurrences 
locate_steps = 10000 #nr of steps after detection
max_ensemble = 50 #nr of initial conditions

#---------------------------system parameters-----------------------------
N = 4 
#N = 3 
g = g_N * (N-1) 
#g = 5.0
#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N
ws = [wc1+0.02:0.01:wc2-0.01;] 
w = wc1 + (wc2-wc1)/2 
w = 0.55
diffeq = (alg = Rodas5(),dt = 1e-3,adaptive=true,reltol = 1e-8, abstol = 1e-10,maxiters=1e15)
#diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,maxiters=1e15)
ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)

#---------------------find two attractors to show as examples----------------------------

#test trajectory for efficient phase space grid choice
traj,_ = trajectory(ds,T_test;Ttr=Ttr,Δt=1e-3)
plot(traj[:,1],traj[:,2])
#need max finding along each dim
max_values = maximum(Matrix(traj);dims=1) #ugly max finding

grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps)

u1 = generate_initial_condition(N)
u2 = generate_initial_condition(N)
u2[1] = 2.0; u2[3] = 1.0
mapper(u1)
mapper(u2)

for i in 1:max_ensemble
    @show mapper(rand(6))
end

attractors = extract_attractors(mapper)
plot_attractors(attractors,access=[1,3,5])
att_plots = []


pl1 = plot()
plot!(pl1,attractors[3][:,1],attractors[3][:,3])
pl2 = plot()
plot!(pl2,attractors[1][:,1],attractors[1][:,3],attractors[1][:,5])
plot(pl1,pl2,layout=(1,2))

for i in 1:length(attractors)
    plot!(pl2,attractors[i][:,1],attractors[i][:,3],attractors[i][:,5])
    #push!(att_plots,pl)
end
plot!(pl2,camera = (100, 20))
#=
nr_attractors = 0
i=0
while (nr_attractors < factorial(N-1) ) || (i < max_ensemble)
    attr_ID = mapper(rand(2N) .* 10)
    nr_attractors = length(mapper.bsn_nfo.attractors)
    @show nr_attractors
    @show i
    @show attr_ID
    i+=1
end
=#