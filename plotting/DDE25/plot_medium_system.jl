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



function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end

N = 10
g = 30.0
wc2 = (N-1)/N 
wc1 = (N-1)/g 
w = 0.4

p = [g,w,N]

u0 = generate_initial_condition(N)
#u0 = rand(2N)
ds = create_ds(u0,p)

traj,t = trajectory(ds,200;Ttr=0,Δt=0.01)
timeseries = Matrix(traj[:,[1:2:2N;]])
pl = plot(t,timeseries[:,1:end],legend=false,xlims=(100,110),la=0.4,ylims=(0,12))

logocolors = Colors.JULIA_LOGO_COLORS
plot!(t,timeseries[:,1:3],lw=3,lc=[logocolors.red logocolors.green logocolors.purple])
plot!(xlabel=L"t",ylabel=L"x_i",guidefontsize=36,
    tickfontsize=20,
    size=(1000,500),left_margin=6Plots.mm,bottom_margin=6Plots.mm,
    dpi=300)

savefig(pl,figs_dir*"system_N$N"*"_w$w"*"_g$g"*".png")