using Plots
using LaTeXStrings
using DelimitedFiles

include("plot_N10_escape_lifetimes.jl")
include("plot_ensemble_average_lyapunov.jl")

plot!(pl_lyap,left_margins=12Plots.mm)
annotate!(pl_lyap, 0.17, 0.95, text(L"(a)", :gray10, :center, 25))
annotate!(pl_esc[1], 0.0004, 2*1e5, text(L"(b)", :gray10, :center, 25))


pl = plot(pl_lyap,pl_esc,layout=(1,2),size=(1200,400),
    right_margins=8Plots.mm,
    bottom_margins=8Plots.mm,dpi=300)

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"transient_chaos_N10_g45_ens1e3_Tmax1e6_ensemble_lyap_escape_times_logdiffparam_cutoff.png")