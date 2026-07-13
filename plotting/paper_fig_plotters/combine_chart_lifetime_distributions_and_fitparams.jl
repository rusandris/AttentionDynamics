using Plots,LaTeXStrings
include("defaults.jl")

#combine distribution plots on the r-curve and quadratic fit params
include("combine_chart_lifetime_distributions.jl")
include("plot_lifetime_distributions_curve_B2B4.jl")


pl_distribs_bfit = plot(pl_distribs_N300, pl_bfit_Tx,layout=(1,2),size=(700,300),title="",
    guidefontsize=18,tickfontsize=10,left_margin=5Plots.mm,bottom_margin=5Plots.mm,
    legendfontsize=12)

#add fig annotations
annotate!(pl_distribs_bfit[1],0.03, 0.07, text(L"(a)", :gray10, :right, annot_size))
annotate!(pl_distribs_bfit[2], 3, 18.5, text(L"(b)", :gray10, :right, annot_size))

savefig(pl_distribs_bfit,"../../paper_figs/lifetime_chart_distribs_quad_fitparams_Tx.png")