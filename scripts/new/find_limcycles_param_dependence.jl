using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using Plots,LaTeXStrings
cd(@__DIR__)

#phase space contraction rate
#trace of the jacobian

function contraction_rate(u,p)
	g,w,N = p
	sum_x = sum(u[1:2:end])
	sum_y = sum(u[2:2:end])
	return g*N - g* sum_y - g*w*(sum_x-sum_y) - N
end

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

grid_size = 20 #phase space grid cell size

T_test = 1000
Ttr=10000
cons_recurrences = 500; 
locate_steps = 10000;
max_ensemble = 50

#---------------------------system parameters-----------------------------
N = 4 # N=4 case
#N = 5 
g = g_N * (N-1) 
wc1 = (N-1)/g
wc2 = (N-1)/N
ws = [wc1+0.02:0.01:wc2-0.01;] 
w = wc1 + (wc2-wc1)/2 

ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)
@show g

#------------------------------------attractor plots------------------------------------
g_val = current_parameters(ds)[1]
set_parameter!(ds,2,0.63)
w_val = current_parameters(ds)[2]
traj,ts = trajectory(ds,200;Ttr=10000,Δt=1e-2)

λ = lyapunov(ds,5000,Ttr=10000,Δt=0.1)

pp1 = plot(traj[:,1],traj[:,2],xlabel=L"x_1",ylabel=L"y_1",label=L"\lambda=%$λ",legendfontsize=14,title=L"N=%$N,g=%$g_val"*L", w=%$w_val",dpi=300)
crs = [contraction_rate(traj[i],current_parameters(ds)) for i in 1:length(traj)]
crs_avg = mean(crs)
pl_cr1 = plot(ts[1:1000],crs[1:1000],xlabel=L"t",ylabel=L"tr(J)",label="",lc=:gray10,ylims=(-Inf,15))
hline!(pl_cr1,[crs_avg],ls=:dash,lc=:red,lw=2,label=L"<tr(J)>=%$crs_avg",legendfontsize=14)
#savefig(pp1,"phase_portrait_N7_example1.png")

set_parameter!(ds,2,0.7)
w_val = current_parameters(ds)[2]
traj,ts = trajectory(ds,200;Ttr=10000,Δt=1e-2)
λ = lyapunov(ds,5000,Ttr=10000,Δt=0.1)

crs = [contraction_rate(traj[i],current_parameters(ds)) for i in 1:length(traj)]
crs_avg = mean(crs)
pl_cr2 = plot(ts[1:1000],crs[1:1000],xlabel=L"t",ylabel=L"tr(J)",label="",lc=:gray10,ylims=(-Inf,15),yticks=[-30:5:15;])
hline!(pl_cr2,[crs_avg],ls=:dash,lc=:red,lw=2,label=L"<tr(J)>=%$crs_avg",legendfontsize=14)
pp2 = plot(traj[:,1],traj[:,2],xlabel=L"x_1",ylabel=L"y_1",label=L"\lambda=%$λ",legendfontsize=14,title=L"N=%$N,g=%$g_val"*L", w=%$w_val",dpi=300)
pp = plot(pp1,pp2,pl_cr1,pl_cr2,layout=(2,2),size=(1200,1200),guidefontsize=18,left_margin=8Plots.mm)
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
	while (nr_attractors < factorial(N-1) ) || (i < max_ensemble)
		attr_ID = mapper(rand(2N) .* 10)
		nr_attractors = length(mapper.bsn_nfo.attractors)
		@show nr_attractors
		@show i
		@show attr_ID
		i+=1
	end
	push!(attractor_finding_results,nr_attractors)

end
