using Revise
using DynamicalSystemsBase
using Attractors
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using Plots,LaTeXStrings #,GLMakie
#plotly()
cd(@__DIR__)
include("../plot_functions/orbit_diagrams.jl")
include("defaults.jl")
#include("../../../scripts/new/bif_floq_od_method_parameters.jl")
figs_dir = "./poster_figs/"
mkpath(figs_dir)


#-----------initial conditions-----------------
#define well-controlled initial condition for each system size

function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end



#parameters from phase diagram 
#gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
#gs = gs_N10 ./(10-1)
#g_N = gs[30]


#--------------------parameters for attractor finding via recurrences---------------
grid_size = 500 #phase space grid cell size
T_test = 1000 #test trajectory for grid edges 
Ttr = 1000 #transient time
cons_recurrences = 500 #nr of recurrences 
locate_steps = 10000 #nr of steps after detection
max_ensemble = 50 #nr of initial conditions



#-------------------------N=3 smaller systems-------------------------------
#system parameters
#N = 4 
N = 3 
g = 3.5 #10 #3.5 #5.0 #g_N * (N-1) 
#g = 4.0
#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N
ws = [wc1+0.02:0.01:wc2-0.01;] 
#w = wc1 + (wc2-wc1)/2 
w = 0.6
#diffeq = (alg = Rodas5(),dt = 1e-3,adaptive=true,reltol = 1e-8, abstol = 1e-10,maxiters=1e15)
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,reltol = 1e-8, abstol = 1e-10,maxiters=1e15)
ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)


#---------------------find two attractors to show as examples----------------------------

#test trajectory for efficient phase space grid choice
traj,_ = trajectory(ds,T_test;Ttr=Ttr,Δt=1e-3)
Plots.plot(traj[:,1],traj[:,2])
#need max finding along each dim
max_values = maximum(Matrix(traj);dims=1) #ugly max finding

grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps)

u1 = generate_initial_condition(N)
u2 = generate_initial_condition(N)
u2[1] = 2.0; u2[3] = 0.5
#u1=rand(6)
#u2=rand(6)
mapper(u1)
mapper(u2)

#=
attractors = extract_attractors(mapper)
#plot_attractors(attractors,access=[1,3],framestyle=:auto)
atr1 = attractors[1]
atr2 = attractors[2]
=#

#plot time series
traj1,t1 = trajectory(ds,100,u1;Ttr=100)
traj2,t2 = trajectory(ds,100,u2;Ttr=100)

#---------------------------plot attractors---------------------------
ms=10
pp = Plots.plot(traj1[:,1],traj1[:,3];
    lw=8,legend_columns=2,
    tickfontsize=tick_fontsize_large+5,
    guidefontsize=guide_fontsize_large+5,
    legendfontsize=legend_fontsize,
    lc=:gray10,xlabel=L"x_1",ylabel=L"x_2",
    label="limit cycle " * L"1")
Plots.plot!(pp, traj2[1:2:end, 1], traj2[1:2:end, 3], lw=2, 
    label="limit cycle " * L"2",
    lc=:gray10,ls=:solid,markerstrokewidth=0.0,
    leftmargin=1Plots.mm,
    bottommargin=5Plots.mm)


#plot timeseries as inset 
lws = [6.0 3.0 0.5]
lw = 3
ls = [:solid :dash :solid] 
colors = collect(Colors.JULIA_LOGO_COLORS)[[1, 2, 4]]

plot!(inset = (1, bbox(0.5,0.15,0.5,0.2)), subplot = 2)
plot!(pp[2],t1[1:200],Matrix(traj1[1:200,[1,3,5]]);
    lc=colors',lw=6,
    ms=0,mc=:purple,markerstrokewidth=0.0,ylabel=L"x_i",legend=false,
    xformatter=:none,ls=:solid,
    labels=[L"x_1" L"x_2" L"x_3"],DEFAULTS...,tickfontsize=tick_fontsize,
    guidefontsize=guide_fontsize+2)


annotate!(pp[2],103, 1.5, text(L"x_1", :gray10, :center, 20))
annotate!(pp[2],107, 1.5, text(L"x_2", :gray10, :center, 20))
annotate!(pp[2],111, 1.5, text(L"x_3", :gray10, :center, 20))


plot!(inset = (1, bbox(0.5,0.4,0.5,0.2)), subplot = 2)
plot!(pp[3],t2[1:200],Matrix(traj2[1:200,[1,3,5]]);
    lc=colors', la=1, 
    ls=:solid,lw=2,
    labels=[L"x_1" L"x_2" L"x_3"],legend=false,
    xlabel=L"t",ylabel=L"x_i",DEFAULTS...,tickfontsize=tick_fontsize,
    guidefontsize=guide_fontsize+2)

annotate!(pp[3],103.5, 1.5, text(L"x_3", :gray10, :center, 20))
annotate!(pp[3],107.5, 1.5, text(L"x_2", :gray10, :center, 20))
annotate!(pp[3],111.5, 1.5, text(L"x_1", :gray10, :center, 20))

plot!(size=(800,700))
pp
#savefig(pp,figs_dir*"N3_multistable_inset.png")

#=    
#=
ts1 = Plots.plot(t1[1:1:100],Matrix(traj1[1:1:100,[1,3,5]]);lc=:purple,la=1,markershape=:circle,
    ms=0,mc=:purple,markerstrokewidth=0.0,ylabel=L"x_i",
    ls=:solid,xformatter=:none,lw=[5.0 3.0 1.0],
    labels=[L"x_1" L"x_2" L"x_3"],DEFAULTS...)
#ls=[:solid :dash :dot]
#plot attractor markers
plot!(ts1,t1[1:5:100],Matrix(traj1[1:5:100,[1,3,5]]);st=:scatter,markershape=:circle,
    ma=1,mc=:purple,markerstrokewidth=0.0,label="",ms=[8 4 2])

ts2 = Plots.plot(t2[1:1:100],Matrix(traj2[1:1:100,[1,3,5]]);lc=:green,
    la=1,markershape=:utriangle,
    ms=0,mc=:green,markerstrokewidth=0.0,
    labels=[L"x_1" L"x_2" L"x_3"],lw=[5.0 3.0 1.0],
    xlabel=L"t",ylabel=L"x_i",DEFAULTS...)
=#

#ls=[:solid :dash :dot]
#plot attractor markers
plot!(ts2,t2[1:5:100],Matrix(traj2[1:5:100,[1,3,5]]),st=:scatter,markershape=:utriangle,
    ma=1,mc=:green,markerstrokewidth=0.0,label="",ms=[8 4 2])

pl_ts = Plots.plot(ts1,ts2,layout=(2,1))

plot!(pl_ts,pp,inset = (1, bbox(0.8,0.0,0.4,1.0)), subplot = 2)

pl = Plots.plot(pp,pl_ts,layout=(1,2),size=fig_size,bottom_margin=5Plots.mm)
savefig(pl,figs_dir*"N3_multistable.png")
=#