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
data_dir_lt_root = "../../data/lifetimes_precomp_orig/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "horizontal_and_curve_article_test_k10_Δt0.02_T500/"
#data_dir_lt = data_dir_lt_root * "horizontal_and_curve_article_test_k10_Δt0.1_T5000/"

chart_times = readdlm(data_dir_lt*"chart_times.txt") 
chart_times = chart_times[[1,3,6]]
T = 500
#plot params
par_labels = [L"B1", L"B2", L"B3", L"B4", L"C"][[1,3,4,5,2]]
colors = [:cyan :lime :blue :darkgreen :gray10][[1, 3, 4, 5, 2]]
markershapes = [:rect :circle :utriangle :ltriangle :rtriangle][[1, 3, 4, 5, 2]]
nbins = 30 #[25,10]
plot_fits = true

include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N300 = pl_loglog_all
pl_linlin_N300 = pl_linlin_all
pl_loglin_N300 = pl_loglin_all
plots_linhist_N300 = plots_linhist

#simple histograms
hist_xlims = [(0,50),(0,1000),(0,2000),(0,2000),(0,2000),
    (0, 20), (0, 50), (0, 200), (0, 200), (0, 200),
    (0, 20), (0, 20), (0, 20), (0, 20), (0, 20)]

for i in 1:length(plots_linhist_N300)
 plot!(plots_linhist_N300[i],xlims=hist_xlims[i])
end

pl_histN300 = plot(plots_linhist_N300..., layout=(3, 5), size=(1200, 600), 
    dpi=300, legendfontsize=15, ylabel="count", xlabel=L"\tau",
    bottom_margin=5Plots.mm, left_margin=8Plots.mm,
    ytickfontsize=5)

plot!(plots_loglog[1][1], left_margin=5Plots.mm)
plots_uniform_binning_N300 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]

#-------------------small system (N=10)---------------------------------
@info "Plotting chart lifetime distributions (N=10)"
#basic parameters 
N = 10
k = 4

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../data/lifetimes_precomp_rescaled/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "/" * "wg_horizontal_and_curve_B1B2B3B4C_k4_Δt0.01_T20000/"
chart_times = readdlm(data_dir_lt*"chart_times.txt")
chart_times = chart_times[[2,4,7]]

#plot params
nbins = 30 #[25,10]
plot_fits = false

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
 
#combine N=300, N=10 plots 
pl = plot(plots_uniform_binning_N10..., plots_uniform_binning_N300...,
    layout=(2,3),size=(1200,600),legend=:topright,bottom_margin=8Plots.mm,
    legendfontsize=12, tickfontsize=10, guidefontsize=18, dpi=300)


#----------------------------save plots---------------------------------------

savedir_pl_root = "../../paper_figs/"
savedir_pl = savedir_pl_root 
#mkpath(savedir_pl)
savefig(pl, savedir_pl * "chart_lifetime_distributions_N10_N300_diff_chart_times_Δt=0.02.png")
savefig(pl_histN300, savedir_pl * "lifetime_chart_distribs_linhists_k$(k)_T$(T)_N300_Δt=0.02.png")
savefig(pl_histN10, savedir_pl * "lifetime_chart_distribs_linhists_k$(k)_T$(T)_N10_Δt=0.01.png")