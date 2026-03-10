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

pl = plot(plots_limc[1:2]..., pl_hopf, plots_limc[3],
    plots_limc[4:5]..., pl_clinic, plots_limc[6],
    layout=l, size=(1600, 600), legend_columns=-1,
    tickfontsize=12, guidefontsize=18, bottom_margin=10Plots.mm, 
    legendfontsize=15, dpi=300, left_margin=10Plots.mm)

figs_dir = "../../paper_figs/"
savefig(pl, figs_dir * "limitcycles_distance_to_fixedpoints_scaling_N2_g5.png")