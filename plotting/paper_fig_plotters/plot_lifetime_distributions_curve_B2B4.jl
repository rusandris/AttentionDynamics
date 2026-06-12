using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
using CurveFit

cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lifetimes_precomp_N300/B24_change_r_Δt0.02_T1000_Ttr2000_k10/" 
#order_pars = [2,9,3,4,5,6,7,1,8]
#order_pars = [1:9;]
order_pars = [17, 1, 16:-1:6..., 18, 5:-1:2...] #reverse([1:18;])
data_dirs = readdir(data_dir_lt_root) #[order_pars]
#files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times] #select binary only

#charting time 
Δt = 2.0 #0.2 #0.02 

#markersize
ms = 5
N = 300
k = 10

pl_distribs = plot()
#colors = [:gray10, :gray30, :gray50, :gray70, :gray90]
colors = range(colorant"lime", stop=colorant"purple", length=length(data_dirs))
#par_labels = [L"r = 14.352 = r_{B_3}", L"r = 10 = r_{B_2}", L"r=9", L"r=8", L"r=6", L"r=7", L"r=5", L"r = 4.45 = r_{B_4}", L"r=3.5"]
par_labels = [data_dirs[i][5:end] for i in 1:length(data_dirs)]
ylabel = L"\log(p_{\Delta t}(\tau))"

quadratic(a, x) = @. a[1] + a[2] * x + a[3] * x^2
cutoff = 30
bs = []
#loop through parameters (lifetime dirs)
for (i, dir) in enumerate(data_dirs)

    #read lifetime file
    files = readdir(data_dir_lt_root*dir)
    #@show files
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
    
    #cutoff infs
    cutoff = findfirst(x -> isinf(x), log10.(counts_discrete))
    if !(isnothing(cutoff))
        bins_discrete = bins_discrete[1:cutoff-1]
        counts_discrete = counts_discrete[1:cutoff-1]
    end

    #--------------------plot uniform discrete binning-------------------------- 
    color = colors[i]
    plot!(pl_distribs, log10.(bins_discrete), log10.(counts_discrete),
        markershape=:circle, #st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.5,
        lc=color,
        alpha=1,
        dpi=300,
        title=L"N = %$(N), k = %$(k), \Delta t = %$(Δt)",
        ms=ms,
        label=par_labels[i],
        xlabel=L"\log(\tau)",
        ylabel=ylabel)


    #try to fit 
    #quadratic 
    n_τ = length(bins_discrete)
    interval = 1:round(Int, 0.99 * n_τ)
    param_guesses = [-1.0, 1.0, -1.0]
    prob = NonlinearCurveFitProblem(quadratic, param_guesses, log10.(bins_discrete[interval]), log10.(counts_discrete[interval]))
    sol = solve(prob)
    ps_fit = sol.u

    push!(bs,ps_fit[end])

    @show ps_fit
    xs_fit = log10.(bins_discrete[interval[1]]:0.1:bins_discrete[interval[end]])
    plot!(xs_fit, quadratic(ps_fit, xs_fit),
        ls=:dot, lc=:gray50, lw=2, label="")
end

figs_dir = "../../paper_figs/"

#plot b values 
#TODO: PUT EXACT VALUES HERE
#r_B4 = 4.456
#r_B3 = 14.352
#r_B2 = 10
#rs = [14.352,10.0:-1.0:7.0...,6.5:-0.25:4.75...,4.456,4.25:-0.25:3.5...]

rs = [3.75,4.0,4.25,4.75,5.25,5.5,5.75,6.25,6.5]

pl_bfit = plot(rs, abs.(bs), markershape=:circle,lc=:gray10,ms=3,label="",
    mc=:gray10, xlabel=L"r", ylabel=L"|b|")

pl = plot(pl_distribs, pl_bfit, layout=(1, 2), size=(1200, 500), 
    left_margin=8Plots.mm, right_margin=8Plots.mm, 
    bottom_margin=8Plots.mm, guidefontsize=25,tickfontsize=16)

#savefig(pl,figs_dir*"lifetime_distributions_and quad_params_B2_to_B4_change_r_charting_time_$(Δt).png")