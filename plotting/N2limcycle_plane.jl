using AttentionDynamics
using OrdinaryDiffEq
using Plots
using LaTeXStrings
using LinearAlgebra
import DynamicalSystems: set_parameter!

include("defaults.jl")
include("plot_functions/timeseries_phase_portrait.jl")
x_fp(w,N,M) = 1/(1+w*M/(N-1)-w*N/(N-1))


u0 = rand(4)
p = [5,0.3,2]
g,w,N = p
Ttr = 100
Δt = 1e-2
T = 40
logocolors = Colors.JULIA_LOGO_COLORS

w_values = [0.3,0.495]
Ts = [40,200]
pp_lims = [(0.0,3.0),(0.0,4.0)]
pp_ticks = [[0.0,1.5,3.0],[0.0,2.0,4.0]]
ts_ticks = [[0.0,1.5,3.0],[0.0,2.5,5.0]]
ts_lims = [(0.0,3.0),(0.0,5.0)]
velm_lims = [(1.0,4.0),(0.0,5.5)]
velm_ticks = [[1.0,2.5,4.0],[0.0,2.5,5.5]]


ds = create_ds(u0,p)
diffeq = (alg=Tsit5(),adaptive=true,dtmax=1e-2,dt=1e-3)

for i in 1:length(w_values)
	w = w_values[i]
	T = Ts[i]
	lims=pp_lims[i]
	ticks=round.(pp_ticks[i];digits=1)
	velm_yticks = velm_ticks[i]
	ts_ylims=ts_lims[i]
	velm_ylims=velm_lims[i]
	ts_yticks = ts_ticks[i] 
	set_parameter!(ds,2,w)
	
	#calc cycle
	cycle = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = u0,p=p,diffeq = diffeq)

	#----------------------phase space velocity---------------------------
	velocity = zeros(4)
	velocity_magnitudes = zeros(size(cycle)[1])

	for i in 1:size(cycle)[1]
		ds.f(velocity,cycle[i,:],p,0) #in-place call to calculate derivatives for every state
		velocity_magnitudes[i] = norm(velocity)
	end

	#plot ts
	ts = plot_timeseries(cycle;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=2,ylabel=L"$x_1,x_2$",yticks=ts_yticks,ylims=ts_ylims,xticks=[])
	
	#plot velocity 
	velm_ts = plot_timeseries(velocity_magnitudes;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=2,lc=:gray10,ylabel=L"$\left|\dot{\mathbf{s}}\right|$",xlabel=L"t",ylims=velm_ylims,yticks=velm_yticks)
	
	#combine the two
	ts_comb = plot(ts,velm_ts,layout=(2,1))
	
	#plot phase portrait projection
	pp = plot(cycle[:,1],cycle[:,3];DEFAULTS...,lc=:gray10,xlabel=L"$x_1$",ylabel=L"$x_2$",lw=3,
	lims=lims,la=0.9,ticks=ticks,label="limit cycle")


	#---------------degeneracy plane------------------
	x2(x1) = 2 - x1 
	plot!(pp,[0,2],x2.([0,2]),lc=:gray30,ls=:dot,lw=3,label="degeneracy plane")


	#---------------------fixed points----------------------

	fp_c = [1.0,1.0]
	fp_x1 = [x_fp(w,2,1),0]
	fp_x2 = [0,x_fp(w,2,1)]

	# coexistence fp
	plot!(pp,[fp_c[1]],[fp_c[2]],st=:scatter,mc=:white,markerstrokecolor=:red,markerstrokewidth=2,ma=1,ms=5,label=L"s_c^*",legend=true,legendfontsize=12)

	#dominance fp
	plot!(pp,[fp_x1[1],fp_x2[1]],[fp_x1[2],fp_x2[2]],st=:scatter,mc=:white,markerstrokecolor=:blue,markerstrokewidth=2,ma=1,ms=5,label=L"s_d^*",legend=true,legendfontsize=12)

	savefig(ts_comb,"timeseries_pvelocity_degen_plane_w=$w"*".pdf")
	savefig(pp,"phase_portrait_degen_plane_w=$w"*".pdf")

end

