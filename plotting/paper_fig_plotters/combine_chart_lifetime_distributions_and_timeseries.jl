using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
import FractalDimensions: linreg
using LinearAlgebra
using CurveFit
cd(@__DIR__)

#-------------------------combine chart lifetime distributions---------------------------------

include("combine_chart_lifetime_distributions.jl")
include("plot_timeseries_to_chart_lifetimes.jl")

#plot!(plots_ts_N300[end],bottom_margin=15Plots.mm)
#plot!(plots_ts_N10[end], bottom_margin=15Plots.mm)
plot!(pl_distribs_N10, titlefontsize=25)
plot!(pl_distribs_N300, titlefontsize=25)

#add panel letter annotations 
annotate!(pl_distribs_N10, -0.1, 0.5, text(L"(a)", :gray10, :center, 25))
annotate!(pl_distribs_N300, 0.0, 0.5, text(L"(b)", :gray10, :center, 25))
annotate!(plots_ts_N10[1], 995, 0.19, text(L"(c)", :gray10, :center, 25))
annotate!(plots_ts_N300[1], 1997, 0.055, text(L"(d)", :gray10, :center, 25))

l = @layout [a{0.3h} b{0.3h}; grid(4, 1) grid(4, 1)]
pl_dist_ts = plot(pl_distribs_N10, pl_distribs_N300, plots_ts_N10..., plots_ts_N300...; 
    layout=l, size=(1200, 1000),left_margin=5Plots.mm,legendfontsize=17)

savefig(pl_dist_ts, savedir_pl * "chart_lifetime_distributions_and_timeseries_N10_N300_orig_params_diff_chart_times_Δt=0.02_T1000_Ttr2000_main_plus.png")


#add one more fig only with B2 and B4 (N=300)
pl_B2 = deepcopy(plots_ts_N300[2])
pl_B4 = deepcopy(plots_ts_N300[4])

plot!(pl_B2,yticks=[0.0,0.2]) #edit ticks
plot!(pl_B4,yticks=[0.0,0.1])
pl_ts_B2B4 = plot(plots_ts_N300[2],plots_ts_N300[4],layout=(2,1),size=(1400,400),
    ylabel=L"L_i",left_margin=10Plots.mm,bottom_margin=12Plots.mm,dpi=300)


savefig(pl_ts_B2B4,"../../paper_figs/timeseries_B2B4_illustration.png")

#side-by-side
pl_B2_sbs = deepcopy(plots_ts_N300[2])
pl_B4_sbs = deepcopy(plots_ts_N300[4])

plot!(pl_B2_sbs,yticks=[0.0,0.1,0.2]) #edit ticks
plot!(pl_B4_sbs,yticks=[0.0,0.05,0.1])
plot!(pl_B2_sbs,ylabel=L"L_i",left_margin=10Plots.mm)

pl_ts_B2B4_sbs = plot(pl_B2_sbs,pl_B4_sbs,layout=(1,2),size=(1400,300),
    xlabel=L"t_L",left_margin=12Plots.mm,bottom_margin=14Plots.mm,dpi=300)

savefig(pl_ts_B2B4_sbs,"../../paper_figs/timeseries_B2B4_illustration_sbs.png")