using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lifetimes_precomp_orig/B24_change_r_Δt0.02_T1000_Ttr2000_k10_Δt0.1_T5000/"
data_dirs = reverse(readdir(data_dir_lt_root))[[8,2,3,4,5,6,1,7]]

#files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times] #select binary only

#charting time 
Δt = 2.0 #0.02 

#markersize
ms = 6
N = 300
k = 10

pl = plot()
colors = [:gray10, :gray30, :gray50, :gray70, :gray90]
colors = range(colorant"lime", stop=colorant"purple", length=length(data_dirs))
par_labels = [L"r=10", L"r=9", L"r=8", L"r=6", L"r=6", L"r=5",L"r=r_{B_4}", L"r=3.5"]
ylabel = L"log(p_{\Delta t}(\tau))"
#loop through parameters (lifetime dirs)
for (i, dir) in enumerate(data_dirs)

    #read lifetime file
    files = readdir(data_dir_lt_root*dir)
    @show files
    file_lifetime = files[findfirst(f -> occursin("charting_$(Δt)", f), files)]#[[13,14,15,16,1:12...]] #select binary only
    @show file_lifetime
    lifetimes = reinterpret(Int64, read(data_dir_lt_root * dir * "/" * file_lifetime))


    #--------------------------get binning for lifetimes--------------------------

    #preset bins, one for each possible \tau value (discrete)
    min_lt, max_lt = extrema(lifetimes)
    bin_edges = min_lt:1:max_lt

    #hist
    h_discrete = fit(Histogram, Float64.(lifetimes), bin_edges)
    h_discrete = normalize(h_discrete; mode=:pdf)
    bins_discrete = h_discrete.edges[1][1:end-1]
    counts_discrete = h_discrete.weights

    #--------------------plot uniform discrete binning-------------------------- 
    color = colors[i]
    plot!(pl, log10.(bins_discrete), log10.(counts_discrete),
        markershape=:circle, st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.1,
        lc=color,
        alpha=1,
        dpi=300,
        title=L"N = %$(N), k = %$(k), \Delta t = %$(Δt)",
        ms=ms,
        label=par_labels[i],
        xlabel=L"\log(\tau)",
        ylabel=ylabel)



end
figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"lifetime_distributions_B2_to_B4_change_r_charting_time_2.png")