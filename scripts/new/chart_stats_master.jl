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

#precompute trajectories
savedir_traj_root = "../../data/trajs_precomp_orig/"
include("precompute_trajectories.jl")

#precompute lifetimes
data_dir_traj = savedir_traj
savedir_lt_root = "../../data/lifetimes_precomp_orig/"
include("precompute_lifetimes_rescaled.jl")

#plot and save distributions 
savedir_pl_root = "../../paper_figs/lifetime_distributions/"
data_dir_lt = savedir_lt
include("../../plotting/paper_fig_plotters/plot_lifetime_distributions_integrated_rescaled.jl")
