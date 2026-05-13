using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lyapunovs_effect_size/"
data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_T_WLC1000_T_fp10_Ttr_WLC1000_Ttr_fp10_dtmax_0.001/"

files = readdir(data_dir_lt)
files_lyaps = readdir(data_dir_lt)[1:end-1][[4,5,6,1,2,3]]
Ns = Int.(readdlm(data_dir_lt * "system_sizes.txt"))

#plot all on one fig 

pl = plot(ylabel=L"\lambda", xlabel=L"w",guidefontsize=20,tickfontsize=15,
    legendfontsize = 15)
for (i,file) in enumerate(files_lyaps)
    
    #read lyap exps 
    lyap_data = readdlm(data_dir_lt * file)
    ws = lyap_data[:,1]
    λs =  lyap_data[:,2]
    plot!(pl, ws, λs, lw=2, label=L"N = %$(Ns[i])")

end
pl