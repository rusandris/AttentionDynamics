using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lyapunovs_effect_size/"
data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws461_T_WLC1000_T_fp10_Ttr_WLC10000_Ttr_fp10_dtmax_0.001/"


files = readdir(data_dir_lt)
files_lyaps = readdir(data_dir_lt)[1:end-1][[4,5,6,1,2,3]]
Ns = Int.(readdlm(data_dir_lt * "system_sizes.txt"))

#plot all on one fig 

pl = plot(ylabel=L"\lambda", xlabel=L"w",guidefontsize=20,tickfontsize=15,
    legendfontsize = 12,legend=:bottomleft,legendcolumns=3,ylims=(-1.0,1.5),
    size=(1200, 400), bottom_margin=8Plots.mm, right_margin=8Plots.mm, 
    left_margin=8Plots.mm,dpi=300)

colors = distinguishable_colors(length(Ns)+3)[3:end]

vline!([0.25], ls=:dash, lc=:gray10, label=false)

for (i,file) in enumerate(files_lyaps)
    N = Ns[i]
    #read lyap exps 
    lyap_data = readdlm(data_dir_lt * file)
    ws = lyap_data[:,1]
    λs =  lyap_data[:,2]
    plot!(pl, ws, λs, lc = colors[i], lw=2, label=L"N = %$(N)")
    vline!([(N - 1) / N], ls=:dash, lc=colors[i], label="")

end
savefig(pl,"../../paper_figs/lyapunovs_system_size_effect_chaotic_width_Ttr1e4.png")