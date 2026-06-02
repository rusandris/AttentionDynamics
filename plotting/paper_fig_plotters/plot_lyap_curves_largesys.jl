using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
using DataFrames, CSV
using StatsBase
using CurveFit  
cd(@__DIR__)

data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_N300_logsys_extended/"
data_dirs = readdir(data_dir_root) #[[2,4]]
data_lyaps1 = readdlm(data_dir_root * data_dirs[1] * "/lyaps_N300.txt")
data_lyaps2 = readdlm(data_dir_root * data_dirs[2] * "/lyaps_N300.txt")
#sort data
data_sorted_lyaps1 = data_lyaps1[sortperm(data_lyaps1[:,2]),:]
data_sorted_lyaps2 = data_lyaps2[sortperm(data_lyaps2[:,2]),:]

Λs1 = data_sorted_lyaps1[:,end]
Λs2 = data_sorted_lyaps2[:,end]

ws1 = data_sorted_lyaps1[:,2]
ws2 = data_sorted_lyaps2[:,2]

#read alphas 
data_alphas1 = readdlm(data_dir_root * data_dirs[1] * "/gs_alphas.txt")
data_alphas2 = readdlm(data_dir_root * data_dirs[2] * "/gs_alphas.txt")

αs = [data_alphas1[2],data_alphas2[2]]

pl_lyaps = plot(xlabel=L"w",ylabel=L"\Lambda_{\log} / \alpha",guidefontsize=18,tickfontsize=15,ylims=(0,8.0),
    title=L"N=300",legend=:topright,legendfontsize=12,dpi=300)
plot!(pl_lyaps, ws1, Λs1 ./ αs[1], label=L"\frac{1}{\varepsilon(N-1)} = 4.45",lw=2,lc=:purple)
plot!(pl_lyaps, ws2, Λs2 ./ αs[2], label=L"\frac{1}{\varepsilon(N-1)} = 6.68",lw=2,lc=:lime)

#approximate values for B2 and B4 levels 
w_B2 = 0.264
w_B4 = 0.446
w_B2_idx = 44
w_B4_idx = 93
w_B2_approx = ws2[w_B2_idx]
w_B4_approx = ws1[w_B4_idx]
@assert isapprox(w_B4,w_B4_approx;atol=1e-2)
@assert isapprox(w_B2,w_B2_approx,atol=1e-2)

ms = 12
plot!([w_B4_approx],[Λs1[w_B4_idx]/αs[1]],st=:scatter,ms=ms,mc=:purple,markershape=:rect,label=L"B_4")
plot!([w_B2_approx],[Λs2[w_B2_idx]/αs[2]],st=:scatter,ms=ms,mc=:lime,markershape=:circle,label=L"B_2")

