include("plot_multistability.jl")
include("plot_attractor_finding.jl")
cd(@__DIR__)

figs_dir = "./poster_figs/"
mkpath(figs_dir)

pl = plot(pp,pl_atr,layout=(1,2),size = (1600,700),
    leftmargin=8Plots.mm,
    bottommargin=10Plots.mm)

savefig(pl,figs_dir*"multistability_combined.png")