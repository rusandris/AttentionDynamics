using DelimitedFiles,Plots,LaTeXStrings

include("defaults.jl")
include("plot_functions/orbit_diagrams.jl")


function read_od(filename)
	od_data = readdlm(filename)
	filtered_od_data = [filter(x -> x != "",od_data[i,:]) for i in 1:size(od_data)[1]]
end



od_datafile = "od_data_N3g_5.0_w_0.6_0.64_tf100_Ttr1000vars1.txt"
od_datafilename,ext = splitext(od_datafile)


od_data = read_od("data/"*od_datafile)
od_wvalues = readdlm("data/od_wvalues_N3g_5.0_w_0.6_0.64_tf100_Ttr1000vars1.txt")

xlims=(0.6,0.64)
xticks = [od_wvalues[1],(od_wvalues[1]+od_wvalues[end])/2,od_wvalues[end]]

od = plot_orbit_diagram(od_wvalues,od_data;OD_DEFAULTS...,xlims=xlims,xticks=xticks)

savefig(od,"figs/orbit_diagram_"*od_datafilename*".png")
od

