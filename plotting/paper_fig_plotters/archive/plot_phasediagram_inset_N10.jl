using Plots,AttentionDynamics,DelimitedFiles,LaTeXStrings

include("defaults.jl")
include("plot_functions/lyapunov.jl")


output_filename = "data/phase_diagram_N10_T=500_Ttr=3000_res_40x40.pdf"
Λ = readdlm("data/N=10.0_T=500_Ttr=3000lyaps.txt")
wvalues = vec(readdlm("data/N=10.0_wvalues.txt"))
gvalues = vec(readdlm("data/N=10.0_gvalues.txt"))

gl = length(gvalues)
wl = length(wvalues)
N = 10

#replace_infs!(Λ,-1)

xticks = [wvalues[1],(wvalues[1]+wvalues[end])/2,(N-1)/N,wvalues[end]]
yticks = [round(gvalues[3]/N;digits=2),(gvalues[1]+gvalues[end])/2/N,gvalues[end]/N]

pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...,xticks=xticks,yticks=yticks,colorbar_ticks=[-1/2,0,1/2])

pd
savefig(pd,"figs/phase_diagram_N10_T=500_Ttr=3000_res_40x40.pdf")
