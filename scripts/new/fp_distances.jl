using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using LinearAlgebra
using StatsBase
#using Attractors
using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

#IDEA: calculate the distance to fixed points along the trajectory

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 3
g = 5.0
w = 0.58
p = [g,w,N]
wc1 = (N-1)/g
wc2 = (N-1)/N

#create ds
u0 = rand(2N) 
ds = create_ds(u0,p;diffeq=diffeq)


#------------------------low w value--------------------------
w = wc1 + 0.01
set_parameter!(ds,2,w)

#generate all fixed single-winner points (with M = 1)
x_fp(w,M,N) = (N-1)/(N-1 + w*M - w*N)
x_sw_val = x_fp(w,1,N)
x_sw = [x_sw_val,x_sw_val,zeros(2N-2)...]
fixedpoints_sw = [zeros(2N) for i  in 1:N]

for (i,fp) in enumerate(fixedpoints_sw)
    fp[2*i-1:2*i] .= x_sw_val
end

#calculate trajectory
traj,t = trajectory(ds,15;Ttr=200,Δt=1e-2)

#calculate distances for every fixed point along traj

dists_fps = []
for (i,fp) in enumerate(fixedpoints_sw)
    dists = [norm(traj[i] .- fp) for i in 1:length(traj)]
    push!(dists_fps,dists)
end

#plots
labels = [L"x_1" L"x_2" L"x_3"]
x_vars = Matrix(traj[:][:,1:2:6]) 
pl_ts = plot(t,x_vars, ylims = (0,4),
    xformatter=:none,ylabel=L"x(t)",labels=labels,
    title=L"w=%$w")

fp_labels = [L"u_1^{*(1)}" L"u_1^{*(2)}" L"u_1^{*(3)}"]
pl_d = plot(t,stack(dists_fps),ylims=(0,5),
    xlabel=L"t",ylabel="distance to fp",labels=fp_labels,
    legendfontsize=20,title=L"w=%$w")

#plot!(pl_ts,pl_d,layout=(2,1),size=(400,600),
#    guidefontsize=18,legendfontsize=15,tickfontsize=10)

#------------------------high w value--------------------------

w = wc2 - 0.01
set_parameter!(ds,2,w)

#generate all fixed single-winner points (with M = 1)
x_fp(w,M,N) = (N-1)/(N-1 + w*M - w*N)
x_sw_val = x_fp(w,1,N)
x_sw = [x_sw_val,x_sw_val,zeros(2N-2)...]
fixedpoints_sw = [zeros(2N) for i  in 1:N]

for (i,fp) in enumerate(fixedpoints_sw)
    fp[2*i-1:2*i] .= x_sw_val
end

#calculate trajectory
traj,t = trajectory(ds,55;Ttr=200,Δt=1e-2)

#calculate distances for every fixed point along traj

dists_fps = []
for (i,fp) in enumerate(fixedpoints_sw)
    dists = [norm(traj[i] .- fp) for i in 1:length(traj)]
    push!(dists_fps,dists)
end

#plots
labels = [L"x_1" L"x_2" L"x_3"]
x_vars = Matrix(traj[:][:,1:2:6]) 
pl_ts2 = plot(t,x_vars, ylims = (0,4),
    xformatter=:none,ylabel=L"x(t)",labels=labels,
    title=L"w=%$w")

fp_labels = [L"u_1^{*(1)}" L"u_1^{*(2)}" L"u_1^{*(3)}"]
pl_d2 = plot(t,stack(dists_fps),ylims=(0,8),
    xlabel=L"t",ylabel="distance to fp",labels=fp_labels,
    legendfontsize=20,title=L"w=%$w")

pl = plot(pl_ts,pl_ts2,pl_d,pl_d2,layout=(2,2),size=(800,600),
    guidefontsize=18,legendfontsize=15,tickfontsize=10,dpi=300)

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"distance_to_fixedpoints_N3_g5.png")