using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

#include helpers
include("defaults.jl")
include("../plot_functions/linreg.jl")

#limit cycles N=2 
include("limcycles_N2_supp.jl")

#scaling 
include("scaling_laws_supp.jl")

#combine plots generated from the scripts
#l = @layout [a{0.8w} [grid(2,1)] ]
#l = (2,4)
#l = [a{0.2w} a{0.2w} a{0.2w}; a{0.2w} a{0.2w} a{0.2w}];
#l = @layout [ grid(2,3) a{1.0h,0.3w}]
l = grid(2,4)

#add fig label annotations 

#top row
annotate!(plots_limc[1], 496, 2.1, text(L"(a)", :gray10, :center, 18))
annotate!(plots_limc[2], 497, 3.1, text(L"(c)", :gray10, :center, 18))
annotate!(plots_limc[3], -1, 2.6, text(L"(e)", :gray10, :center, 18))
annotate!(pl_hopf, -0.03, 1.53, text(L"(g)", :gray10, :center, 18))


#bottom row
annotate!(plots_limc[4], 450, 5.2, text(L"(b)", :gray10, :center, 18))
annotate!(plots_limc[5], 458, 8.2, text(L"(d)", :gray10, :center, 18))
annotate!(plots_limc[6], -1.5, 4, text(L"(f)", :gray10, :center, 18))
annotate!(pl_clinic, -3.5, 2.7, text(L"(h)", :gray10, :center, 18))


pl = plot(plots_limc[1:3]..., pl_hopf,
    plots_limc[4:6]..., pl_clinic,
    layout=l, size=(1600, 600), legend_columns=-1,
    tickfontsize=12, guidefontsize=18, bottom_margin=10Plots.mm, 
    legendfontsize=15, dpi=300, left_margin=10Plots.mm)

figs_dir = "../../paper_figs/"
savefig(pl, figs_dir * "limitcycles_distance_to_fixedpoints_scaling_N2_g5.png")