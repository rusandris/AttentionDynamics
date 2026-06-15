using Plots

#combine distribution plots on the r-curve and quadratic fit params
include("plot_lifetime_distributions_curve_B2B4.jl")
include("combine_chart_lifetime_distributions.jl")

pl_distribs_bfit = plot(pl_distribs_N300, pl_bfit,layout=(1,2),size=(700,300),title="",
    guidefontsize=18,tickfontsize=10,left_margin=5Plots.mm,bottom_margin=5Plots.mm,
    legendfontsize=12)

#add fig annotations
annotate!(pl_distribs_bfit[1],0.05, 0.0, text("(a)", :gray10, :right, 15))
annotate!(pl_distribs_bfit[2],2.0, 50.0, text("(b)", :gray10, :right, 15))

savefig(pl_distribs_bfit,"../../paper_figs/lifetime_chart_distribs_quad_fitparams.png")