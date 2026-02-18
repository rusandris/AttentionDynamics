using Plots,AttentionDynamics,DelimitedFiles,LaTeXStrings

include("defaults.jl")
include("plot_functions/lyapunov.jl")


output_filename = ("figs/phase_diagram_N4.0_wvalues_wstart_0.67_wstop_0.71_T1000_Ttr15000_res_40x40.pdf")
Λ = readdlm("data/N4.0_wstart_0.67_wstop_0.71_gstart_7_gstop_12_T1000_Ttr15000lyaps_res_40x40.txt")
wvalues = vec(readdlm("data/N4.0_wvalues_wstart_0.67_wstop_0.71_T1000_Ttr15000_res_40x40.txt"))
gvalues = vec(readdlm("data/N4.0_gvalues_gstart_7_gstop_12_T1000_Ttr15000_res_40x40.txt"))
bp_points = readdlm("data/N4bifpoints.txt")


gl = length(gvalues)
wl = length(wvalues)
N = 4

xticks = [wvalues[1],round((N-1)/N;digits=2),wvalues[end]]
yticks = [round(gvalues[3]/N;digits=2),round((gvalues[1]/N+gvalues[end]/N)/2;digits=2),round(gvalues[end]/N;digits=2)]



replace_infs!(Λ,-1)
pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...,xticks=xticks,yticks=yticks)


for bp in eachrow(bp_points)
	bptype = Symbol(bp[2])
	plot!([bp[3]],[bp[1] ./ N],st=:scatter,ms=1,mc=BIF_COLORS[bptype])
end


savefig(pd,output_filename)
