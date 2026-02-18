using Plots,LaTeXStrings
using DelimitedFiles

include("plot_functions/bifurcation_diagrams.jl")
include("defaults.jl")

fp2_bifpoints = readdlm("data/fp2_bifpoints.txt")

pl = plot_bifdiagram("data/fp2_contdata.txt","data/fp2_bifpoints.txt";BIF_DEFAULTS...,xlabel=L"w",ylabel=L"x_1^{*}",ylims=(-0.01,5),xlims=(-0.01,1.01),xticks=[0.0,0.5,1.0],lc=:gray10)
plot_bifdiagram!(pl,"data/po_hopfN2_bifdata.txt",:Stable;BIF_DEFAULTS...,lc=:red)

plot_bifdiagram!(pl,"data/fp1_contdata.txt",:Unstable;BIF_DEFAULTS...)
plot_bifdiagram!(pl,"data/fp3_contdata.txt","data/fp3_bifpoints.txt";BIF_DEFAULTS...,lc=:gray10)
#plot_bifdiagram!(pl,"data/fp4_contdata.txt","data/fp4_bifpoints.txt";BIF_DEFAULTS...)


# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(0.5,5,0.5,0), opacity=.1,c=:blue)
plot!(rectangle(fp2_bifpoints[2],5,0,0), opacity=.1,c=:red)
annotate!(0.1,4.5,text("coexistence", 12,:right,rotation=90))
annotate!(0.35,4.5,text("periodic", 12,:right,rotation=90))
annotate!(0.9,4.5,text("dominance",12,:right,rotation=90))


annotate!(0.8,1.3,text(L"s^*_{c}", 18,:gray10,:center))
annotate!(0.8,3.8,text(L"s^*_{d}", 18,:gray10,:center))
plot!(title=L"N=2")



savefig(pl,"figs/bifdiagN2_g4.pdf")
