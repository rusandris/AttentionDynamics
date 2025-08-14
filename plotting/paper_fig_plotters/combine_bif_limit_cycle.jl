#combine bifurcation plot and limit cycle plot 

include("N2limcycle_plane.jl")
include("plot_small_system_bif.jl")

figs_dir = "../../paper_figs/"
mkpath(figs_dir)

l = @layout [a{0.7w} b{0.3w} ]
pl_comb = plot(pl_bif,pl_limc,layout = l,size=(1300,600),left_margin=5Plots.mm,bottom_margin=5Plots.mm)
savefig(pl_comb,figs_dir*"pl_N2_bif_limc_combined_test.pdf")