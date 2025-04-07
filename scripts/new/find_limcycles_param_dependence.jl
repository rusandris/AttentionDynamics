using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

grid_size = 200 #phase space grid cell size

T_test = 1000
Ttr=10000
cons_recurrences = 500; 
locate_steps = 10000;
max_ensemble = 50
attractor_finding_results = []

#---------------------------system parameters-----------------------------
#N = 4 # N=4 case
N = 5 
g = g_N * (N-1) 
wc1 = (N-1)/g
wc2 = (N-1)/N
ws = [wc1+0.02:0.01:wc2-0.01;] 
w = wc1 + (wc2-wc1)/2 

ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)
@show g

#attractor plots
g_val = current_parameters(ds)[1]
set_parameter!(ds,2,ws[1])
w_val = current_parameters(ds)[2]
traj,_ = trajectory(ds,200;Ttr=10000,Δt=1e-2)
pp1 = plot(traj[:,1],traj[:,2],xlabel=L"x_1",ylabel=L"y_1",title=L"N=%$N,g=%$g_val"*L", w=%$w_val",dpi=300)
#savefig(pp1,"phase_portrait_N7_example1.png")

set_parameter!(ds,2,ws[10])
w_val = current_parameters(ds)[2]
traj,_ = trajectory(ds,200;Ttr=10000,Δt=1e-2)
pp2 = plot(traj[:,1],traj[:,2],xlabel=L"x_1",ylabel=L"y_1",title=L"N=%$N,g=%$g_val"*L", w=%$w_val",dpi=300)
pp = plot(pp1,pp2,layout=(1,2),size=(1200,600),guidefontsize=18,left_margin=8Plots.mm)
savefig(pp,"../../phase_portrait_N$(N)_examples.png")

attractor_finding_results = []

for (i,w) in enumerate(ws)
	@show w
	set_parameter!(ds,2,w)
	
	#test trajectory for efficient phase space grid choice
	traj,_ = trajectory(ds,T_test;Ttr=Ttr)
	#need max finding along each dim
	max_values = maximum(Matrix(traj);dims=1) #ugly max finding

	grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])
	
	mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
	consecutive_recurrences = cons_recurrences, attractor_locate_steps = locate_steps)
	nr_attractors = 0
	i=0
	while (nr_attractors < factorial(N-1) ) || (i > max_ensemble)
		attr_ID = mapper(rand(2N) .* 10)
		nr_attractors = length(mapper.bsn_nfo.attractors)
		@show attr_ID
		i+=1
	end
	push!(attractor_finding_results,nr_attractors)

end
