using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
import FractalDimensions: linreg
using LinearAlgebra
using CurveFit
cd(@__DIR__)

#-------------------------combine chart lifetime distributions---------------------------------


#-------------------MAIN--big system (N=300)---------------------------------
@info "Plotting chart lifetime distributions (N=300)"
#basic parameters 
N = 300
k = 10

#lifetime paths 
#data_dir_lt = savedir_lt
data_dir_lt_root = "../../data/lifetimes_precomp_N300/"
dirs_lifetimes = readdir(data_dir_lt_root)
data_dir_lt = data_dir_lt_root * "B2B3B4_large_charting_k10_Δt0.02_T1000_Ttr2000/"
#data_dir_lt = data_dir_lt_root * "horizontal_and_curve_article_test_k10_Δt0.02_T500/"

chart_times = readdlm(data_dir_lt*"chart_times.txt") 
T = 1000

#--------------------------plot params---------------------------
logocolors = Colors.JULIA_LOGO_COLORS
#plot_order = [1, 2, 4, 5, 3]
#plot_order = [1, 3, 4, 2]
#plot_order = [1:4;]
plot_order = [1:3;]
fileset_order = [1,2,3] #[2, 3, 4] #B2,B3,B4
#fileset_order = [3,5,6,7,1,4]

#par_labels = [L"B_1", L"B_2", L"B_3", L"B_4"][plot_order]
par_labels = [L"B_2", L"B_3", L"B_4"][fileset_order]
par_labels_string = ["B2", "B3", "B4"][fileset_order]

#colors = [:cyan :gray60 :lime :orange :purple][plot_order]
#colors = [:cyan :lime :orange :purple][plot_order]
colors = [:lime :orange :purple][fileset_order]
#markershapes = [:diamond, :star4 ,:circle, :utriangle, :rect][plot_order]
#markershapes = [:diamond, :circle, :utriangle, :rect][plot_order]
markershapes = [:dtriangle, :utriangle, :rect][fileset_order]
nbins = 30 #[25,10]

#fit details, setting 
plot_fits = true
linfit_idx = 1
quadfit_idxs = [1,2,3]
chart_time_fit_idx = [1,2]
fit_lengths = [1.0,1.0]
msw=0.5

include("../../plotting/paper_fig_plotters/plot_lifetime_distributions.jl")

pl_loglog_N300 = pl_loglog_all
pl_linlin_N300 = pl_linlin_all
pl_loglin_N300 = pl_loglin_all
plots_linhist_N300 = plots_linhist

plot!(plots_loglog[1][1], left_margin=5Plots.mm)
plots_uniform_binning_N300 = [plots_loglog[i][1] for i in 1:length(plots_loglog)]

#---------------------------select 2 plots from many (main)--------------------------------------

N = 300
pl_distribs_N300 = plot(plots_uniform_binning_N300...,
    layout=(1,2),legend=:topright, ylims=(-4.5,0.1),
    bottom_margin=10Plots.mm, left_margin=10Plots.mm,
    ylabel=L"\log(p(\tau))", 
    xlims = [0.2,1.7],xticks = [0.5,1.0,1.5],
    tickfontsize=15, guidefontsize=22, dpi=300,
    size=(1200,400),legendfontsize=15)


#test 
