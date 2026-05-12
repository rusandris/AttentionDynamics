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
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c_article/"
params_orig_horizontal = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled_horizontal = readdlm(param_dir * "params_w_g.txt")
params_changed_horizontal = readdlm(param_dir * "params_c.txt")
params_changed_horizontal = vec(params_changed_horizontal)
cs = vec(params_changed_horizontal)


#alpha (rp=rc) direction on phasediagram
param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_alpha_N.txt")
K, c, α0, N_orig = params_orig_curve
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = vec(params_changed_curve)


#set up small system (N=10) with the original parameters 
#so it matches the rescaled phase diagram

N_orig = params_orig_horizontal[end]
N_small = 10
ws_horizontal = params_rescaled_horizontal[:, 1]
gs_horizontal = params_rescaled_horizontal[:, 2] * (N_small - 1)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2] * (N_small - 1)

αs_horizontal = rs[2] ./ ((1 ./ ws_horizontal .- 1) .* (cs .* K .* (N_small - 1))) #all the same
α = αs_horizontal[1] #constant
αs_curve = rs ./ gs_curve #all the same

#gather g,w pairs (horizontal,curve)
p_sets_r = [[rs[i],rs[i]] for i in 1:length(rs)]
p_sets_c = [[cs[i]] for i in 1:length(cs)]
p_sets_c = [[3.4]]


#--------------------------------precompute trajectories-------------------------------------
#settings of precomputed trajs
N = N_small
T = 5000 #10000 #100 #5000
Ttr = 5000 #500
Δt_traj = 1e-2 #5*1e-2
nr_ics = 10 #number of initial conditions
ic_counter = 0

#phasediag_ids = ["B1", "B2", "C", "B3", "B4"]
phasediag_ids = ["C"]

#paths
savedir_traj_root = "../../data/trajs_precomp_orig_N10/"
savedir_traj = savedir_traj_root * "B1B2B3B4C_Δt$(Δt_traj)_T$(T)_Ttr$(Ttr)/"
mkpath(savedir_traj)

include("precompute_trajectories.jl")


#-------------------------------precompute lifetimes------------------------------------------
#lifetime settings 
k = 10
chart_times = [0.02,0.1,0.2,0.5,1.0,2.0,3.0,5.0] #[2.0]
#chart_times = [0.02]

#path
data_dir_traj = "../../data/trajs_precomp_orig_N10/B1B2B3B4C_Δt0.01_T5000_Ttr5000/"
savedir_lt_root = "../../data/lifetimes_precomp_orig_N10/"
savedir_lt = savedir_lt_root * "B1B2B3B4Cmod34_N10_Δt0.01_T5000_Ttr5000/"
mkpath(savedir_lt)
writedlm(savedir_lt * "chart_times.txt", chart_times)

#par_names = ["B1", "B2", "C", "B3", "B4"]
par_names = ["C"]
include("precompute_lifetimes.jl")


