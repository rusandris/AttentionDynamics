using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
using CurveFit
include("defaults.jl")

cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lifetimes_precomp_N300/lifetimes_rcurve_allvalues_B24_change_r_Δt0.02_T1000_Ttr2000_k10/"
#order_pars = [2,9,3,4,5,6,7,1,8]
#order_pars = [1:9;]
order_pars = [17, 1, 16:-1:6..., 18, 5:-1:2...] #reverse([1:18;])
data_dirs = readdir(data_dir_lt_root)[order_pars]
#files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times] #select binary only

rs = [14.352, 10.0:-1.0:7.0..., 6.5:-0.25:4.75..., 4.456, 4.25:-0.25:3.5...]
#rs = [3.75,4.0,4.25,4.75,5.25,5.5,5.75,6.25,6.5]

#charting time 
Δt = 2.0 #0.2 #0.02 

#-----------------------setup plot-----------------------

#markersize
ms_default = 2
ms_high = 8 
N = 300
k = 10

ylabel = L"\log(p_{\Delta t}(\tau))"

pl_distribs = plot(dpi=300, left_margin=8Plots.mm,
    xlabel=L"\log(\tau)",ylabel=ylabel) #title=L"N = %$(N), k = %$(k), \Delta t = %$(Δt)"
#colors = [:gray10, :gray30, :gray50, :gray70, :gray90]
colors = range(colorant"lime", stop=colorant"purple", length=length(data_dirs)-1)
pushfirst!(colors,colorant"orange")
colors = reverse(colors)
#par_labels = [L"r = 14.352 = r_{B_3}", L"r = 10 = r_{B_2}", L"r=9", L"r=8", L"r=6", L"r=7", L"r=5", L"r = 4.45 = r_{B_4}", L"r=3.5"]
#par_labels = [data_dirs[i][5:end] for i in 1:length(data_dirs)]
par_labels = [L"r = r_{B3} = %$(round(rs[1];digits=2))", 
    L"r = r_{B2} = %$(rs[2])", 
    L"r = r_{B4} = %$(round(rs[end-4];digits=2))"]

markershape_default =:circle

#fitting function
quadratic(a, x) = @. a[1] + a[2] * x + a[3] * x^2
cutoff = 30
bs = []

#loop through parameters (lifetime dirs)
for (i, dir) in enumerate(reverse(data_dirs))

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
    if i == length(data_dirs)
        par_label = par_labels[1]  
        markershape =:utriangle
        ms = ms_high
    elseif i == length(data_dirs) - 1
        par_label = par_labels[2]
        markershape =:dtriangle
        ms = ms_high
    elseif i == 5
        par_label = par_labels[3]
        markershape =:rect
        ms = ms_high
    else 
        par_label = ""
        markershape = markershape_default
        ms = ms_default
    end
    
    plot!(pl_distribs, log10.(bins_discrete), log10.(counts_discrete),
        markershape=markershape, #st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.5,
        lc=color,
        alpha=1,
        ms=ms,
        label=par_label)


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
        ls=:dot, lc=color, lw=2, label="")
end

figs_dir = "../../paper_figs/"

#plot b values 
#TODO: PUT EXACT VALUES HERE
#r_B4 = 4.456
#r_B3 = 14.352
#r_B2 = 10

#plot quadratic coefficients in terms of r
bs = reverse(bs)
pl_bfit_r = plot(rs[1:end-1], abs.(bs[1:end-1]), markershape=:circle, lc=:gray10, ms=3, label="",
    mc=:gray10, xlabel=L"r", ylabel=L"|b|",ylims=(-5,50)) #yaxis=:log10, xaxis=:log10)


#overlay markers of B2, B3, B4 
plot_B3 = true
plot!([rs[2]], [abs.(bs[2])], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label=L"B_2")
plot_B3 && plot!([rs[1]], [abs.(bs[1])], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label=L"B_3")
plot!([rs[end-4]], [abs.(bs[end-4])], st=:scatter,markershape = :rect, mc=:purple, ms=8, label=L"B_4")


#plot quadratic coefficients in terms of w
c0 = 2.4;α0 = 0.005
wg_pairs = rescaled_params.(rs,rs,1.0,c0,α0,N)
ws = [wg_pairs[i][1] for i in 1:length(wg_pairs)]
pl_bfit_w = plot(ws[1:end-1], abs.(bs[1:end-1]), markershape=:circle, lc=:gray10, ms=3, 
    mc=:gray10, xlabel=L"w", ylabel=L"|b|", ylims=(-5,50),label="") #yaxis=:log10, xaxis=:log10)

##overlay markers of B2, B3, B4 
plot_B3 = false
plot!([ws[2]], [abs.(bs[2])], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label=L"B_2")
plot_B3 && plot!([ws[1]], [abs.(bs[1])], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label=L"B_3")
plot!([ws[end-4]], [abs.(bs[end-4])], st=:scatter,markershape = :rect, mc=:purple, ms=8, label=L"B_4")

#--------------plot supplimentary figure--------------------


#add colorbar 
cbar_range = reverse(rs)
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [3:2:15;], 
    c=colors, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-6Plots.mm,
    tickfontsize=tick_fontsize_large-5,
    right_margin=18Plots.mm)

#combine heatmap with the separate colorbar
l = @layout [a{0.95w} b]
pl_distribs_cbar = plot!(pl_distribs, pl_cbar, layout=l, tick_direction=:out)

pl_distfit_sup = plot(pl_distribs_cbar, pl_bfit_r, layout=(1, 2), size=(1000, 400), 
    legendfontsize=12, bottom_margin=8Plots.mm, top_margin=12Plots.mm, guidefontsize=20,tickfontsize=12)

#add fig annotations

annotate!(pl_distfit_sup[1], 0.05, 0.8, text(L"(a)", :gray10, :center, annot_size))
annotate!(pl_distfit_sup[3], 0.12, 56.0, text(L"(b)", :gray10, :center, annot_size))

savefig(pl_distfit_sup,figs_dir*"lifetime_distributions_and quad_params_B2_to_B4_change_r_charting_time_$(Δt).png")