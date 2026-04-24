using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

#setup system 

#N = 300
N = 10
#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds
x0 = rand(N) * 0.01
u0 = zeros(2N)
u0[1:2:end] .= x0

r = 10.0;
K = 1.0;
c = 2.4;
α0 = 0.005;
p_orig = [r, r, K, c, α0, N]

w_B2 = 0.2640565204592287
g_B2 = 60.20066889632107

w_B2_round = 0.264
g_B2_round = 60.2
#------------------------version 1: (LY system) keep r anc c the same, scale α----------------------------
α = r ./ ((1 ./ w_B2 .- 1) .* (c .* K .* (N - 1))) #all the same
α_round = round(α,digits=5) 
p = [r, r, K, c, α, N]

ds = create_original_ds(u0, p; diffeq=diffeq)
T = 500
traj, t = trajectory(ds, T; Ttr=0, Δt=1e-2)
pl_L = plot(t,Matrix(traj)[:,1:2:end],legend=false,
    title=L"N=%$(N),r=%$(r),K=%$(K),c=%$(c),\alpha \approx %$(α_round)", 
    ylabel=L"L_i", xlabel=L"t", guidefontsize=20)

#--------------------------version 2: (xy system) rescale trajectory itself----------------------------
p_resc = [g_B2,w_B2,N]
η_L, η_Y = rescaling_factors(r, r, K, c, α, N)

u0_resc = zeros(2N)
u0_resc[1:2:end] .= u0[1:2:end] ./ η_L #rescale Ls
u0_resc[2:2:end] .= u0[2:2:end] ./ η_Y #rescale Ys

ds_resc = create_ds(u0_resc, p_resc; diffeq=diffeq)
traj_resc, t_resc = trajectory(ds_resc, T*α, u0_resc; Ttr=0, Δt=1e-2)
pl_resc = plot(t_resc, Matrix(traj_resc)[:,1:2:end],legend=false,
    title=L"N=%$(N),1/\epsilon \approx %$(g_B2_round),w \approx %$(w_B2_round))",
    ylabel=L"x_i",xlabel=L"t^'", guidefontsize=20)

traj_resc_L = Matrix(traj_resc)[:,1:2:end]
traj_resc_L .*= η_L
pl_resc_resc = plot(t_resc ./ α, traj_resc_L,legend=false,
    title=L"N=%$(N),\approx %$(g_B2_round),w \approx %$(w_B2_round))",
    ylabel=L"x_i \cdot \eta_{L}", xlabel=L"t^' / \alpha", guidefontsize=20)

pl = plot(pl_L, pl_resc, pl_resc_resc, layout=(3,1),size=(500,800),left_margin=5Plots.mm,dpi=300)
savefig(pl,"../../paper_figs/check_timescaling_rescaling_N10_B2_T500.png")