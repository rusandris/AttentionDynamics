using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)


#-----------------------------------time params---------------------------------------

#length of precomputed trajs 
T = 100000 #10000 #100 #5000
Ttr = 500 #500
Δt = 1e-1 #5*1e-2

#-----------------------------------setup system--------------------------------------

diffeq2 = (alg=Tsit5(), dt=1e-2, dtmax=1e-1, adaptive=true, maxiters=1e15)

N = 10
g = 18.0
w = 0.7
u0 = rand(2N) * 0.1

ds = create_ds(u0, [g, w, N]; diffeq=diffeq2) #no noise


savedir_root = "../../data/trajs_precomp_rescaled_triangle/"
mkpath(savedir_root)

#-------------------change rp---------------------
#parameters to change 
gs = [16.499999999999996, 20.0]
ws = [0.574468085106383, 0.574468085106383]

p_sets = [[gs[i], ws[i]] for i in 1:length(gs)]

#precompute trajectories
savedir = savedir_root * "rp_ic10/"
mkpath(savedir)
writedlm(savedir * "g_params.txt", gs)
writedlm(savedir * "w_params.txt", ws)

nr_ics = 10

precompute_trajectories(ds, nr_ics, p_sets, [1, 2]; Δt=1e-1,
    T_traj=T, Ttr=Ttr, savedir=savedir, filename="_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat")

#-------------------change alpha---------------------
#parameters to change 
gs = [20.0, 16.499999999999996]
ws = [0.574468085106383, 0.6206896551724138]

p_sets = [[gs[i],ws[i]] for i in 1:length(gs)]

#precompute trajectories
savedir = savedir_root*"alpha_ic10/"
mkpath(savedir)
writedlm(savedir*"g_params.txt", gs)
writedlm(savedir*"w_params.txt", ws)

nr_ics = 10

precompute_trajectories(ds, nr_ics, p_sets, [1,2]; Δt=1e-1,
    T_traj=T, Ttr=Ttr, savedir=savedir, filename="_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat")

#-------------------change competition---------------------
#parameters to change 
gs = [16.499999999999996, 16.499999999999996]
ws = [0.6206896551724138, 0.5744680851063829]

p_sets = [[gs[i], ws[i]] for i in 1:length(gs)]

#precompute trajectories
savedir = savedir_root * "c_ic10/"
mkpath(savedir)
writedlm(savedir * "g_params.txt", gs)
writedlm(savedir * "w_params.txt", ws)

nr_ics = 10

precompute_trajectories(ds, nr_ics, p_sets, [1, 2]; Δt=1e-1,
    T_traj=T, Ttr=Ttr, savedir=savedir, filename="_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat")

#-------------------change r---------------------
#parameters to change 
gs = [14.210526315789473, 18.421052631578945]
ws = [0.6551724137931035, 0.5944380069524914]

p_sets = [[gs[i], ws[i]] for i in 1:length(gs)]

#precompute trajectories
savedir = savedir_root * "rp=rc_ic10/"
mkpath(savedir)
writedlm(savedir * "g_params.txt", gs)
writedlm(savedir * "w_params.txt", ws)

nr_ics = 10

precompute_trajectories(ds, nr_ics, p_sets, [1, 2]; Δt=1e-1,
    T_traj=T, Ttr=Ttr, savedir=savedir, filename="_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat")


