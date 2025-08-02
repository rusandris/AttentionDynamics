using Plots,LaTeXStrings
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using Attractors
cd(@__DIR__)
figs_dir = "../../paper_figs/"

#traj
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.3
p = [g,w,N]
u0 = rand(2N) .* 40
ds = create_ds(u0,p;diffeq=diffeq)

T_plot_chaos = 1000
Ttr_chaos = 4000
Δt = 0.001

#plot transient
traj_chaos,t_tr = trajectory(ds,T_plot_chaos,u0;Δt = Δt)

x1 = traj_chaos[:,1]
x2 = traj_chaos[:,3]
y1 = traj_chaos[:,2]
y2 = traj_chaos[:,4]

plot!(x1 .- x2, y1 .- y2,
    lw=1,la=0.2,
    label="",ylims=(-4,4),xlims=(-25,25))

#IDEA: compute unstable limit cycles for this parameter?
#using Bifkit.jl? 