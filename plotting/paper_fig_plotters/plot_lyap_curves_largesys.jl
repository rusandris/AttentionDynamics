using Plots, LaTeXStrings
using AttentionDynamics
using DelimitedFiles
using DataFrames, CSV
using StatsBase
using CurveFit  
cd(@__DIR__)

data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_N300_logsys_extended/"

data_dirs = readdir(data_dir_root) #[[2,4]]
data_lyaps1 = readdlm(data_dir_root * data_dirs[1] * "/lyaps_N300.txt")
data_lyaps2 = readdlm(data_dir_root * data_dirs[3] * "/lyaps_N300.txt")

#sort data
data_sorted_lyaps1 = data_lyaps1[sortperm(data_lyaps1[:,2]),:]
data_sorted_lyaps2 = data_lyaps2[sortperm(data_lyaps2[:,2]),:]

Λs1 = data_sorted_lyaps1[:,end]
Λs2 = data_sorted_lyaps2[:,end]

ws1 = data_sorted_lyaps1[:,2]
ws2 = data_sorted_lyaps2[:,2]

#read alphas 
data_alphas1 = readdlm(data_dir_root * data_dirs[1] * "/gs_alphas.txt")
data_alphas2 = readdlm(data_dir_root * data_dirs[3] * "/gs_alphas.txt")

αs = [data_alphas1[2],data_alphas2[2]]
ylabel = L"\Lambda" #L"\Lambda_{\log} / \alpha"
label = L"1/(T_x(N-1))"
pl_lyaps = plot(xlabel=L"w",ylabel=ylabel,guidefontsize=18,tickfontsize=15,ylims=(0,8.0),
    legend=:topright,legendfontsize=12,dpi=300,xlims=(0.0,0.9),xticks=[0.1:0.2:0.9;])
plot!(pl_lyaps, ws1, Λs1 ./ αs[1], label = label *  L" = 2.98",lw=2,lc=:purple)
plot!(pl_lyaps, ws2, Λs2 ./ αs[2], label= label *  L" = 6.68",lw=2,lc=:lime)

#approximate values for B2 and B4 levels 
w_B2 = 0.264
w_B4 = 0.446
w_B2_idx = 44
w_B4_idx = 37
w_B2_approx = ws2[w_B2_idx]
w_B4_approx = ws1[w_B4_idx]
@assert isapprox(w_B4,w_B4_approx;atol=1e-2)
@assert isapprox(w_B2,w_B2_approx,atol=1e-2)

ms = 12
plot!([w_B4_approx],[Λs1[w_B4_idx]/αs[1]],st=:scatter,ms=ms,mc=:purple,markershape=:rect,label="")
plot!([w_B2_approx],[Λs2[w_B2_idx]/αs[2]],st=:scatter,ms=ms,mc=:lime,markershape=:circle,label="")

#add annotations
annotate!(w_B2 + 0.05, Λs2[w_B2_idx]/αs[2] + 0.5, text(L"B_2", :lime, :center, 20))
annotate!(w_B4 + 0.08, Λs1[w_B4_idx]/αs[1] + 0.2, text(L"B_4", :purple, :center, 20))


#plot in function of r 
data_dir_root_r = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_N300_logsys_change_r/"
data_lyaps_r = readdlm(data_dir_root_r * readdir(data_dir_root_r)[1] * "/lyaps_N300.txt")
data_sorted_lyaps_r = data_lyaps_r[sortperm(data_lyaps_r[:,2]),:]
r_values = data_sorted_lyaps_r[:,2]

#read const params for plot in function of r
const_params = readdlm(data_dir_root_r * readdir(data_dir_root_r)[1] * "/const_params.txt")
K,c0,α0 = const_params
Λs_r = data_sorted_lyaps_r[:,end] ./ α0


wg_vals = rescaled_params.(r_values,r_values,K,c0,α0,300)
ws_r = [wg_vals[i][1] for i in 1:length(wg_vals)]

#add double axis

plot(r_values, Λs_r, lc=:blue,mc=:blue,markershape=:circle, 
    foreground_color_axis=:blue,foreground_color_border=:blue,
    ylabel=L"\Lambda", xlabel=L"r",guidefontsize=25,ylims=(2,9),
    label="",plot_title=L"N=300,c=2.4,α=0.005",
    plot_titlevspan=0.1,dpi=300)
plot!(twiny(), ws_r, Λs_r, lc=:red, ylims=(2,9), mc=:red, 
    xflip=true, markershape=:circle, xlabel=L"w", guidefontsize=25,
        foreground_color_axis=:red,foreground_color_border=:red,
        label="")

#=
plot(r_values, Λs_r, lc=:blue,mc=:blue,markershape=:circle, 
    foreground_color_axis=:blue,foreground_color_border=:blue,
    ylabel=L"\Lambda", xlabel=L"r",guidefontsize=25,ylims=(2,9),
    label="",plot_titlevspan=0.1,dpi=300,size=(600,400))

plot(ws_r, Λs_r, lc=:red, ylims=(2,9), mc=:red, 
    xflip=true,  ylabel=L"\Lambda", markershape=:circle, xlabel=L"w", guidefontsize=25,
        foreground_color_axis=:red,foreground_color_border=:red,
        label="",size=(600,400),dpi=300)
=#