using AttentionDynamics,DelimitedFiles,Plots,LaTeXStrings

include("defaults.jl")
include("plot_functions/orbit_diagrams.jl")


function read_od(filename)
	od_data = readdlm(filename)
	filtered_od_data = [filter(x -> x != "",od_data[i,:]) for i in 1:size(od_data)[1]]
end



od_datafile = "od_data_N4g_6.0_w_0.69_0.72_tf1000_Ttr1000_ic5vars1.txt"
od_datafilename,ext = splitext(od_datafile)


od_data = read_od("data/"*od_datafile)
od_wvalues = readdlm("data/od_wvalues_N4g_6.0_w_0.69_0.72_tf1000_Ttr1000_ic5vars1.txt")

xlims=(0.699,0.72)
starred_values = [0.7,0.71]
xticks = [starred_values[1],starred_values[2],od_wvalues[end]]

od = plot_orbit_diagram(od_wvalues,od_data;OD_DEFAULTS...,xlims=xlims,ylims=(1.6,3.8),xticks=xticks)
scatter!(starred_values,[3.7,3.7],markershape=:star5,c=[:gray70, :gray40],ms=10,markerstrokecolor=[:gray70, :gray40])


savefig(od,"figs/orbit_diagram_"*od_datafilename*".pdf")
od

