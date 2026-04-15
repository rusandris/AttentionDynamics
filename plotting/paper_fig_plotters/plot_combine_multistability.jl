include("plot_multistability.jl")
include("plot_attractor_finding.jl")
cd(@__DIR__)

figs_dir = "../../paper_figs/"
mkpath(figs_dir)

annotate!(pp, -0.25, 2.8, text(L"(a)", :gray10, :center, 30))
annotate!(pl_atr, -5, 35, text(L"(b)", :gray10, :center, 30))

pl = plot(pp,pl_atr,layout=(1,2),size = (1600,700),
    leftmargin=8Plots.mm,
    bottommargin=10Plots.mm)

savefig(pl,figs_dir*"multistability_combined.png")