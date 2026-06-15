using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
import FractalDimensions: linreg
using LinearAlgebra
using CurveFit
cd(@__DIR__)

#-------------------------combine chart lifetime distributions---------------------------------


#-------------------big system (N=300)---------------------------------
@info "Plotting chart lifetime distributions (N=300)"
#basic parameters 
N = 300
k = 10

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../data/lifetimes_precomp_N300/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "chart_lifetimes_B1B2B3B4C_Cmod34_k10_Δt0.02_T1000_Ttr2000/"
#data_dir_lt = data_dir_lt_root * "horizontal_and_curve_article_test_k10_Δt0.02_T500/"

chart_times = readdlm(data_dir_lt*"chart_times.txt") 
chart_times = chart_times[[1,3,6]]
T = 1000

#--------------------------plot params---------------------------
logocolors = Colors.JULIA_LOGO_COLORS
#plot_order = [1, 2, 4, 5, 3]
#plot_order = [1, 3, 4, 2]
#plot_order = [1:4;]
plot_order = [1:3;]
fileset_order = [2, 3, 4]
#fileset_order = [3,5,6,7,1,4]

#par_labels = [L"B_1", L"B_2", L"B_3", L"B_4"][plot_order]
par_labels = [L"B_1", L"B_2", L"B_3", L"B_4", L"B_{12}", L"B_{24}"][fileset_order]
#colors = [:cyan :gray60 :lime :orange :purple][plot_order]
#colors = [:cyan :lime :orange :purple][plot_order]
colors = [:cyan :lime :orange :purple :gray20 :gray60][fileset_order]
#markershapes = [:diamond, :star4 ,:circle, :utriangle, :rect][plot_order]
#markershapes = [:diamond, :circle, :utriangle, :rect][plot_order]
markershapes = [:diamond, :circle, :utriangle, :rect, :circle, :circle][fileset_order]
nbins = 30 #[25,10]

#fit details, setting 
plot_fits = true
linfit_idx = 1
quadfit_idx = 3


include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N300 = pl_loglog_all
pl_linlin_N300 = pl_linlin_all
pl_loglin_N300 = pl_loglin_all
plots_linhist_N300 = plots_linhist

plot!(plots_loglog[1][1], left_margin=5Plots.mm)
plots_uniform_binning_N300 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]

#=
#-------------------small system (N=10)---------------------------------
@info "Plotting chart lifetime distributions (N=10)"
#basic parameters 
N = 10
k = 4

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../data/lifetimes_precomp_orig_N10/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "/" * "B1B2B3B4C_N10_Δt0.01_T5000_Ttr1000/"
chart_times = readdlm(data_dir_lt*"chart_times.txt")
chart_times = chart_times[[1,3,6]]

#plot params
nbins = 30 #[25,10]
plot_fits = false
fileset_order = [1:4;]

include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N10 = pl_loglog_all
pl_linlin_N10 = pl_linlin_all
pl_loglin_N10 = pl_loglin_all
plots_linhist_N10 = plots_linhist

#simple histograms
hist_xlims = [(0, 100), (0, 300), (0, 300), (0, 400), (20, 200),
    (0, 20), (0, 50), (0, 50), (0, 50), (0, 20),
    (0, 10), (0, 10), (0, 10), (0, 10), (0, 10)]

for i in 1:length(plots_linhist_N10)
    plot!(plots_linhist_N10[i], xlims=hist_xlims[i])
end

pl_histN10 = plot(plots_linhist_N10..., layout=(3, 5), size=(1200, 600), 
    dpi=300, legendfontsize=15,ylabel="count",xlabel=L"\tau",
    bottom_margin=5Plots.mm,
    ytickfontsize=5)

plot!(plots_loglog[1][1],left_margin=5Plots.mm)
plots_uniform_binning_N10 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]
 
#---------------------------combine N=300, N=10 plots-------------------------------------
pl_sup = plot(plots_uniform_binning_N10..., plots_uniform_binning_N300...,
    layout=(2,3),size=(1200,600),legend=:topright,bottom_margin=8Plots.mm,
    legendfontsize=12, tickfontsize=10, guidefontsize=18, dpi=300)

=#
#---------------------------select two figs (main)--------------------------------------

N = 300
pl_distribs_N300 = plot!(plots_uniform_binning_N300[end],legend=:topright, ylims=(-4.5,0.1),
    bottom_margin=2Plots.mm, ylabel=L"\log(p_{\Delta t}(\tau))",  xlims = [0.2,1.7],xticks = [0.5,1.0,1.5],
    tickfontsize=15, guidefontsize=22, dpi=300)

#=
#N = 10
plot!(plots_uniform_binning_N10[end],ylabel=L"\log(p_{\Delta t} (\tau))")

pl_distribs_N10 = plot!(plots_uniform_binning_N10[end],legend=:topright, ylims=(-4.5,0.1),
    bottom_margin=2Plots.mm, xlims = [0.2,1.7],xticks = [0.5,1.0,1.5],
    legendfontsize=12, tickfontsize=15, guidefontsize=22, dpi=300)

pl_main = plot(pl_distribs_N10, pl_distribs_N300,
    layout=(1, 2), size=(600, 300), legend=:topright, ylims=(-4.5,0.1),
    bottom_margin=2Plots.mm, xlims = [0.2,1.7],xticks = [0.5,1.0,1.5],
    legendfontsize=15,ms=5, tickfontsize=10, guidefontsize=18, dpi=300)

=#
#----------------------------save plots---------------------------------------

#savedir_pl_root = "../../paper_figs/"
#savedir_pl = savedir_pl_root 
#mkpath(savedir_pl)
#savefig(pl_main, savedir_pl * "chart_lifetime_distributions_N10_N300_orig_params_diff_chart_times_Δt=0.02_T1000_Ttr2000_main.png")
#savefig(pl_sup, savedir_pl * "chart_lifetime_distributions_N10_N300_orig_params_diff_chart_times_Δt=0.02_T1000_Ttr2000_noC.png")
#savefig(pl_histN300, savedir_pl * "lifetime_chart_distribs_linhists_orig_params_k$(k)_T$(T)_N300_Δt=0.02_T1000_Ttr2000_C34.png")
#savefig(pl_histN10, savedir_pl * "lifetime_chart_distribs_linhists_orig_params_k$(k)_T$(T)_N10_Δt=0.02_T1000_Ttr2000_C34.png")