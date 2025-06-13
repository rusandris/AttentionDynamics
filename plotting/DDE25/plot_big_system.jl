using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../defaults.jl")

figs_dir = "./poster_figs/"
mkpath(figs_dir)

#-------------------------------plot big system--------------------------
diffeq = (alg = Tsit5(),)
N = 300

#try original equations
u0 = generate_initial_condition(N)
rc=9.0;rp=9.0;K=1.0;c=2.4;α=0.005
p = [rc,rp,K,c,α]
ds = create_original_ds(u0,p)

wc1 = (N-1)/gs[1]
wc2 = (N-1)/N

traj,t = trajectory(ds,100;Ttr=1000,Δt=0.01)
timeseries = Matrix(traj[:,[1:2:2N;]])

linecolors = fill(:gray50,N-1)
#linecolors[1] = :gray10
pl = plot(t,timeseries[:,1:end],legend=false,ylims=[0.0,0.02],xlims=(1050,1075),la=0.3)
plot!(t,timeseries[:,50:90],legend=false,ylims=[0.0,0.02],la=1.0,lw=3)
plot!(xlabel=L"t",ylabel=L"L_i",guidefontsize=26,
    xticks=[1055:5:1070;],yticks=[0.0,0.01,0.02],tickfontsize=22,
    size=(1000,500),left_margin=6Plots.mm,bottom_margin=6Plots.mm,
    dpi=300)
savefig(pl,figs_dir*"big_system_N$N"*"_r$rc"*"_K$K"*"_c$c"*"_alfa_$α"*".png")