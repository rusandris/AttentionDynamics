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
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c/"
params_orig = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled = readdlm(param_dir * "params_w_g.txt")
params_changed = readdlm(param_dir * "params_c.txt")

#--------------------------------precompute trajectories-------------------------------------
#settings of precomputed trajs
T = 1000 #10000 #100 #5000
Ttr = 500 #500
Δt_traj = 1e-2 #5*1e-2
nr_ics = 20 #number of initial conditions

#paths
savedir_traj_root = "../../data/trajs_precomp_orig/"
savedir_traj = savedir_traj_root * "c_Δt$(Δt_traj)_T$(T)_Ttr$(Ttr)/"
mkpath(savedir_traj)

include("precompute_trajectories.jl")

#-------------------------------precompute lifetimes------------------------------------------
#lifetime settings 
k = 10
chart_times = [0.01,0.1, 0.5, 1.0] #2.0]

#path
data_dir_traj = savedir_traj
savedir_lt_root = "../../data/lifetimes_precomp_orig/"
savedir_lt = savedir_lt_root * "c_k$(k)_Δt$(Δt_traj)_T$(T)/" #* "c2/" #"rp_rc/"
mkpath(savedir_lt)

include("precompute_lifetimes.jl")

#---------------------------plot and save distributions---------------------------------------
savedir_pl_root = "../../paper_figs/lifetime_distributions/"
data_dir_lt = savedir_lt
include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")
