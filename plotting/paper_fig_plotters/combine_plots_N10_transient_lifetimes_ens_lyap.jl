using Plots
using LaTeXStrings
using DelimitedFiles

include("plot_N10_escape_lifetimes.jl")
include("plot_ensemble_average_lyapunov.jl")

pl = plot(pl_lyap,pl_esc,layout=(1,2),size=(1000,400),
    left_margins=8Plots.mm, right_margins=8Plots.mm,
    bottom_margins=8Plots.mm,dpi=300)

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"transient_chaos_N10_ensemble_lyap_escape_times.png")