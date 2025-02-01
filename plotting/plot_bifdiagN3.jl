using Plots,LaTeXStrings
using DelimitedFiles

include("plot_functions/bifurcation_diagrams.jl")
include("defaults.jl")

fp2_bifpoints = readdlm("data/N3fp2_bifpoints.txt")

#plot_bifdiagram!(pl,"data/po_hopfN2_bifdata.txt",:Stable;BIF_DEFAULTS...)
pl = plot(;BIF_DEFAULTS...,xlabel=L"w",ylabel=L"x_1^{*}",ylims=(-0.1,6),xlims=(0.0,1))
plot_bifdiagram!(pl,"data/N3fp1_contdata.txt",:Unstable;BIF_DEFAULTS...,lc=:gray10)
plot_bifdiagram!(pl,"data/N3fp2_contdata.txt","data/N3fp2_bifpoints.txt";BIF_DEFAULTS...)
plot_bifdiagram!(pl,"data/N3fp3_contdata.txt","data/N3fp3_bifpoints.txt";BIF_DEFAULTS...,lc=:gray10)
plot_bifdiagram!(pl,"data/N3fp4_contdata.txt","data/N3fp4_bifpoints.txt";BIF_DEFAULTS...,lc=:gray10)
plot_bifdiagram!(pl,"data/N3po1_contdata.txt","data/N3po1_bifpoints.txt";BIF_DEFAULTS...,lc=:red)
plot_bifdiagram!(pl,"data/N3po2_contdata.txt","data/N3po2_bifpoints.txt";BIF_DEFAULTS...,lc=:red)
plot_bifdiagram!(pl,"data/N3po_hopf_contdata.txt",:Unstable;BIF_DEFAULTS...,lc=:red)

#plot_bifdiagram!(pl,"data/N3po1_contdata_avg.txt","data/N3po1_bifpoints.txt";BIF_DEFAULTS...)
#plot_bifdiagram!(pl,"data/N3po2_contdata_avg.txt","data/N3po2_bifpoints.txt";BIF_DEFAULTS...,lc=:gray50)

#lens!([0.39,2/3],[0.99,1.02];BIF_DEFAULTS...,inset = (1, bbox(0.1, 0.1, 0.5, 0.2)),lw=1,lc=:gray70,la=0.5)
#plot!(pl[2],ticks=false)


# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(2/3,6,2/3,0), opacity=.1,c=:blue)
plot!(rectangle(fp2_bifpoints[3],6,0,0), opacity=.1,c=:red)
annotate!(0.1,4.5,text("coexistence", 12,:right,rotation=90))
annotate!(0.5,3.5,text("periodic", 12,:center,rotation=90))
annotate!(0.9,4.5,text("dominance",12,:right,rotation=90))

annotate!(0.9,1.35,text(L"s^*_{c}", 18,:gray10,:center))
annotate!(0.8,3.8,text(L"s^*_{d}", 18,:gray10,:center))
plot!(title=L"N=3")

savefig(pl,"figs/bifdiagN3_g5.pdf")
