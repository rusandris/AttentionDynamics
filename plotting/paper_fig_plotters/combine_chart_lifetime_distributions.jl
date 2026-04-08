using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#-------------------------combine chart lifetime distributions---------------------------------


#-------------------big system (N=300)---------------------------------
#basic parameters 
N = 300
k = 10

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../../../../Dropbox/att_dyn/data/data_lifetimes/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * dirs_lifetimes[1] * "/"
#chart_times = readdlm(data_dir_lt*"chart_times.txt") 
chart_times = [0.1,0.5,1.0]

#plot params

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

par_labels = [L"B1", L"B2", L"B3", L"B4", L"C"]

colors = [:cyan :lime :gray10 :blue :darkgreen]
markershapes = [:rect :circle :utriangle :ltriangle :rtriangle]
nbins = 30 #[25,10]

include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N300 = pl_loglog_all
pl_linlin_N300 = pl_linlin_all
pl_loglin_N300 = pl_loglin_all

plot!(plots_loglog[1][1], left_margin=5Plots.mm)
plots_uniform_binning_N300 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]

#-------------------small system (N=10)---------------------------------
#basic parameters 
N = 10
k = 4

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../data/lifetimes_precomp_rescaled/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "/" * "wg_horizontal_and_curve_B1B2B3B4C_k4_Δt0.01_T20000/"

chart_times = readdlm(data_dir_lt*"chart_times.txt")
chart_times = chart_times[[1,3,4]]
#plot params

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

par_labels = [L"B1", L"B2", L"B3", L"B4", L"C"]

colors = [:cyan :lime :gray10 :blue :darkgreen]
markershapes = [:rect :circle :utriangle :ltriangle :rtriangle]

nbins = 30 #[25,10]
include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N10 = pl_loglog_all
pl_linlin_N10 = pl_linlin_all
pl_loglin_N10 = pl_loglin_all

plot!(plots_loglog[1][1],left_margin=5Plots.mm)
plots_uniform_binning_N10 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]
 
#combine N=300, N=10 plots 
pl = plot(plots_uniform_binning_N10..., plots_uniform_binning_N300...,
    layout=(2,3),size=(1200,600),legend=:bottomleft,bottom_margin=8Plots.mm,
    legendfontsize = 12,tickfontsize=10,guidefontsize=18,dpi=300)


#----------------------------save plots---------------------------------------

savedir_pl_root = "../../paper_figs/"
savedir_pl = savedir_pl_root 
#mkpath(savedir_pl)
savefig(pl,savedir_pl * "chart_lifetime_distributions_N10_N300_diff_chart_times.png")