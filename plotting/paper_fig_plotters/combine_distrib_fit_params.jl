include("plot_lifetime_dsitributions_curve_B2B4.jl")
include("combine_chart_lifetime_distributions.jl")

pl = plot(pl_distribs_N300, pl_bfit,layout=(1,2),size=(800,400))