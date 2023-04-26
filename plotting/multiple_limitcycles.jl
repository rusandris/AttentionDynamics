using AttentionDynamics
using OrdinaryDiffEq
using Plots
using LaTeXStrings

include("defaults.jl")
include("plot_functions/timeseries_phase_portrait.jl")

u0 = rand(6)
p = [4,0.55,3]
g,w,N = p
Ttr = 100
Δt = 1e-2
T = 50
logocolors = Colors.JULIA_LOGO_COLORS
legendlabels = [L"$x_1$" L"$x_2$" L"$x_3$"]
pp_labels = (x = L"$x_1$",y= L"$x_2$")
ts_labels = (x = L"$t$",y= L"$x_i$")

ds = create_ds(u0,p)
diffeq = (alg=Tsit5(),adaptive=false,dt=1e-3)

state1 = [0.810765070284594, 0.7733829745894634, 0.6229954171573494, 0.5856302648429531, 0.3348947466616211, 0.7743679638708577]
state2 = [0.2849241076093, 0.08481930770740131, 0.4901061113231271, 0.4619221556280978, 0.27540205357693637, 0.4713032328251774]

cycle1 = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = state1,p=p,diffeq = diffeq)
cycle2= evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = state2,p=p,diffeq = diffeq)

#----------------------------plotting----------------------------------

ts1 = plot_timeseries(cycle1;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=4,
	xlims = (Ttr,Ttr+10),ylims=(0,3.5),yticks = [0:1:3.5;],xticks=[],
	ylabel=ts_labels.y,xlabel="",legend=:bottomright,label=legendlabels)

@show state1

ts2 = plot_timeseries(cycle2;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=4,
	xlims = (Ttr,Ttr+10),ylims=(0,3.5),yticks = [0:1:3.5;],xticks=Ttr:3:Ttr+T,
	ylabel=ts_labels.y,xlabel=ts_labels.x)
@show state2


pp = plot_phase_portrait(cycle1[:,1],cycle1[:,3];DEFAULTS...,lc=:gray50,xlabel=pp_labels.x,ylabel=pp_labels.y,xlims = (0,3.5),ylims = (0,3.5),lw=3,xticks=[0,1,2,3],yticks=[0,1,2,3])
plot_phase_portrait!(pp,cycle2[:,1],cycle2[:,3];DEFAULTS...,lc=:gray10,xlabel=pp_labels.x,ylabel=pp_labels.y,xlims = (0,3.5),ylims = (0,3.5),lw=3)

l = @layout [a{0.5w,0.58h} [b{0.5h};c{0.5h}]]
pl = plot(pp,ts1,ts2,layout=l)

savefig(pl,"figs/multiple_limcycles_N3_w0.55_g4.pdf")

