using AttentionDynamics
using Revise
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

#trajectories, lifetimes and distribution plots in one script

#---------------------read parameter values------------------

#horizontal direction on phasediagram 
param_dir_root = "../../../data/parameters/"
param_dir = param_dir_root * "c_article/"
params_orig_horizontal = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled_horizontal = readdlm(param_dir * "params_w_g.txt")
params_changed_horizontal = readdlm(param_dir * "params_c.txt")
params_changed_horizontal = vec(params_changed_horizontal)
cs = vec(params_changed_horizontal)


#alpha (rp=rc) direction on phasediagram
param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_alpha_N.txt")
K, c, α, N_orig = params_orig_curve
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = vec(params_changed_curve)

#equal distance to coexistence curve
#on g/(N-1) graph
#between big and small systems 

N_orig = params_orig_horizontal[end]
N_small = 10
ws_horizontal = params_rescaled_horizontal[:, 1]
gs_horizontal = params_rescaled_horizontal[:, 2] * (N_orig - 1)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2] * (N_orig - 1)

#gather g,w pairs (horizontal,curve)
p_sets_r = [[rs[i], rs[i]] for i in 1:length(rs)]
p_sets_c = [[cs[i]] for i in 1:length(cs)]

#--------------------------------precompute trajectories-------------------------------------
#settings of precomputed trajs
T = 1000 #10000 #100 #5000
Ttr = 2000 #500
Δt_traj = 2 * 1e-2 #5*1e-2
nr_ics = 100 #number of initial conditions
ic_counter = 0

#phasediag_ids = ["B1", "B2", "C", "B3", "B4"]
phasediag_ids = ["B2"]

#paths
savedir_traj_root = "../../../data/trajs_precomp_orig/B1B2B3B4C_RKNONA_Δt$(Δt_traj)_T$(T)_Ttr$(Ttr)/"
savedir_traj = savedir_traj_root * "B2/"
mkpath(savedir_traj)

#setup system 

N = 300

#------------------------------solver settings----------------------------------
dtmax = 1e-3 #5*1e-2
#diffeq = (alg=Tsit5(), dt=1e-2, dtmax=dtmax, adaptive=true, maxiters=1e15)
diffeq = (alg=RK4(), dt=1e-3, dtmax=dtmax, adaptive=false, maxiters=1e15)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds
x0 = rand(N) * 0.01
u0 = zeros(2N)
u0[1:2:end] .= x0

p_init = [rs[2], rs[2], K, c, α, N]
ds = create_original_ds(u0, deepcopy(p_init); diffeq=diffeq)

#write traj specific params
writedlm(savedir_traj * "traj_spec_params.txt", [N Δt_traj T]) #save traj specific params

#select only c value that corresponds to B2
p_B2 = p_sets_c[2:2]

#precompute trajectories (c)
set_parameter!(ds, 1, rs[2])
set_parameter!(ds, 2, rs[2])
precompute_trajectories(ds, nr_ics, p_B2, [4]; Δt=Δt_traj,
    T_traj=T, Ttr=Ttr, savedir=savedir_traj, p_names=["c"],
    ic_counter=ic_counter, phasediag_ids=phasediag_ids,
    filename="_N$(N)_rp$(rs[2])_rp$(rs[2])_K$(K)_alpha$(α)_T$(T)_Ttr$(Ttr)_deltat$(Δt_traj).dat")

