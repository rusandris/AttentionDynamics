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

#alpha (rp=rc) direction on phasediagram
param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_N.txt")
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
#params_changed_curve = readdlm(param_dir * "params_c.txt")

#equal distance to coexistence curve
#on g/(N-1) graph
#between big and small systems 

N_orig = params_orig_horizontal[end]
N_small = 10
ws_horizontal = params_rescaled_horizontal[:, 1]
gs_horizontal = params_rescaled_horizontal[:, 2]  * (N_small - 1)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2] * (N_small - 1)

#gather g,w pairs (horizontal,curve)
p_sets = [[gs_horizontal[i], ws_horizontal[i]] for i in 2:2]
#push!(p_sets, [gs_curve[1], ws_curve[1]], [gs_curve[3], ws_curve[3]])


#--------------------------------precompute trajectories-------------------------------------
#settings of precomputed trajs
T = 20000 #10000 #100 #5000
Ttr = 500 #500
Δt_traj = 1e-2 #5*1e-2
nr_ics = 10 #number of initial conditions
#=
#paths

savedir_traj_root = "../../data/trajs_precomp_rescaled/"
savedir_traj = savedir_traj_root * "wg_horizontal_and_curve_B1B2B3B4C_Δt$(Δt_traj)_T$(T)_Ttr$(Ttr)/"
mkpath(savedir_traj)

#phasediag_ids = ["B2"]
#phasediag_ids = ["B1","B2","B3","C","B4"]
#phasediag_ids = ["C", "B3", "B4"]
include("precompute_trajectories_rescaled.jl")

=#
#-------------------------------precompute lifetimes------------------------------------------
#lifetime settings 
k = 6
chart_times = [0.01, 0.02, 0.1, 0.2, 0.5, 1.0, 2.0] #2.0]
#chart_times = [0.5]

#path
#data_dir_traj = savedir_traj
datadir_traj_root = "../../data/trajs_precomp_rescaled/"
data_dir_traj = datadir_traj_root * "wg_horizontal_and_curve_B1B2B3B4C_Δt0.01_T20000_Ttr500/"
savedir_lt_root = "../../data/lifetimes_precomp_rescaled/"
savedir_lt = savedir_lt_root * "wg_horizontal_and_curve_B1B2B3B4C_k$(k)_Δt$(Δt_traj)_T$(T)/" #* "c2/" #"rp_rc/"
mkpath(savedir_lt)

#p_values = [[1, 2, 3, 4]] #[[3]]#
#p_names = ["B", "C"]#["g", "w"] #["c"] #["rp", "rc"]

writedlm(savedir_lt * "chart_times.txt", chart_times)
#chart_times = [2.0]

include("precompute_lifetimes.jl")

#---------------------------plot and save distributions---------------------------------------
#=
#par labels 
p_names_ltx = [L"B"] #["g", "w"] #[L"r_p", L"r_c"]

#label names  
par_labels = []
for j in 1:length(p_values[1])
    #put strings together
    par_names = ""
    for i in 1:length(p_values)
        par_names *= p_names_ltx[i] * L"_%$(p_values[i][j])"
    end
    @show par_names
    push!(par_labels, par_names)
end
=#
par_labels = [L"B1", L"B2", L"C", L"B3", L"B4"]

colors = [:cyan :lime :gray10 :blue :darkgreen]
markershapes = [:rect :circle :utriangle :ltriangle :rtriangle]

nbins = 30 #[25,10]

#set paths
savedir_pl_root = "../../paper_figs/lifetime_distributions_rescaled/"
savedir_pl = savedir_pl_root * "wg_horizontal_and_curve_article_params_k$(k)_Δt$(Δt_traj)_T$(T)/"
mkpath(savedir_pl)
data_dir_lt = savedir_lt

include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")
