using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
using CurveFit
using DataFrames, CSV
using LsqFit
cd(@__DIR__)
include("defaults.jl")

R_squared(resid,y) = 1 - sum(resid .^ 2) / sum( ( y .- mean(y) ) .^2)    
R_squared_check(sol,y) = 1 - rss(sol;weighted=false)/ (var(y;corrected=false)*length(y))

function fit_irls(model, x, y, p0;abstol=1e-2,maxiter=100)
    n = length(y)
    w_prev = ones(length(y))
    w = deepcopy(w_prev)
    p_prev = deepcopy(p0)
    p = deepcopy(p0)
    #iterate until params, weights have converged
    i = 0
    while i < maxiter 
        i += 1

        fit = curve_fit(model, x, y, w, p)
        p = fit.param

        #reweight
        w_prev .= w
        R = fit.resid .^2
        rss_mean = mean(R)
        w .= [1/(R[i] + rss_mean) for i in 1:n] 
        w .= w ./ sum(w) #normalize

        if norm(p .- p_prev) < abstol && norm(w .- w_prev) < abstol
            return p,w
        end 

        p_prev .= p
    end
    @warn "Did not converge after $maxiter iterations!"
    return p,w
end


figs_dir = "../../paper_figs/"

#path to lifetimes
data_dir_lt_root = "../../data/lifetimes_precomp_N300/lifetimes_rcurve_allvalues_B24_change_r_Δt0.02_T1000_Ttr2000_k10/"
#data_dir_lt_root = "../../data/lifetimes_precomp_N300/lifetimes_rcurve_allvalues_B24_change_r_Δt0.02_T1000_Ttr2000_k10_large_charting_times"

order_of_pars = [17, 1, 16:-1:6..., 18, 5:-1:2...] #reverse([1:18;])
data_dirs_all = reverse(readdir(data_dir_lt_root)[order_of_pars])

#leave out some values off the end (distributions too short)
data_dirs = data_dirs_all[5:end] #[7:end]

#charting time 
chart_times = readdlm(data_dir_lt_root*"/"*data_dirs[1]*"/chart_times.txt") 
#Δt_chart = chart_times[1]
Δt_chart = chart_times[end] 
T = 1000

#hardcode r values
rs_all = [14.352, 10.0:-1.0:7.0..., 6.5:-0.25:4.75..., 4.456, 4.25:-0.25:3.5...]
rs = reverse(rs_all[1:end-4])#reverse(rs_all[1:end-6])

#-----------------------setup plot-----------------------

#markersize
ms_default = 2
ms_high = 8 
N = 300
k = 10

ylabel = L"\log(p(\tau))"

pl_distribs = plot(dpi=300, left_margin=8Plots.mm,
    xlabel=L"\log(\tau)",ylabel=ylabel);

pl_distribs_part = deepcopy(pl_distribs); 
pl_distribs_lin = deepcopy(pl_distribs);

#title=L"N = %$(N), k = %$(k), \Delta t = %$(Δt)"
#colors = [:gray10, :gray30, :gray50, :gray70, :gray90]

pl_distribs_exp = plot(dpi=300, left_margin=8Plots.mm,
    xlabel=L"\log(\tau)",ylabel=ylabel);

colors = range(colorant"purple", stop=colorant"lime", length=length(data_dirs)-1);
push!(colors,colorant"orange");
#par_labels = [L"r = 14.352 = r_{B_3}", L"r = 10 = r_{B_2}", L"r=9", L"r=8", L"r=6", L"r=7", L"r=5", L"r = 4.45 = r_{B_4}", L"r=3.5"]
#par_labels = [data_dirs[i][5:end] for i in 1:length(data_dirs)]
par_labels = [L"r = r_{B3} = %$(round(rs[end];digits=2))", 
    L"r = r_{B2} = %$(rs[end-1])", 
    L"r = r_{B4} = %$(round(rs[1];digits=2))"]

markershape_default =:circle

#fitting function
#quadratic(a, x) = @. a[1] + a[2] * x + a[3] * x^2 #CurveFit
quadratic(x,p) = @. p[1] + p[2] * x + p[3] * x^2 #LsqFit

#fit settings 
param_guesses = [-1.0, 1.0, -1.0]
abstol = 1e-3
maxiter = 300

fitparams_df = DataFrame(parameter=String[],chart_time=Float64[],nr_bins=Int64[],
    min_lt = Int64[], max_lt = Int64[], tot_lts = Int64[], q1=Float64[],
    q2=Float64[],q_ratio=Float64[])

distr_df = DataFrame()
df_len = 50

bins = []
bin_counts = []
R2_vals = []
resid_ss_vals = []
sols = []
R2_lin_vals = []
resid_ss_lin_vals = []
nbins_exp = 8

#loop through parameters (lifetime dirs)
for (i, dir) in enumerate(data_dirs)

    #read lifetime file
    files = readdir(data_dir_lt_root * "/" * dir)
    @show files
    #search for file with the given charting time 
    file_lifetime = files[findfirst(f -> occursin("charting_$(Δt_chart)", f), files)]#[[13,14,15,16,1:12...]] #select binary only
    println("File: " * file_lifetime[1:25] * "...")
    lifetimes = reinterpret(Int64, read(data_dir_lt_root * "/" * dir * "/" * file_lifetime))

    #--------------------------get binning for lifetimes--------------------------

    #preset bins, one for each possible \tau value (discrete)
    min_lt, max_lt = extrema(lifetimes)
    nr_tot_lts = length(lifetimes)
    bin_edges = min_lt:1:max_lt

    #hist
    h_discrete = fit(Histogram, Float64.(lifetimes), bin_edges)
    h_discrete = normalize(h_discrete; mode=:none) #:pdf
    bins_discrete = h_discrete.edges[1][1:end-1]
    counts_discrete = h_discrete.weights
    
    #cutoff infs
    cutoff_inf = findfirst(x -> isinf(x), log10.(counts_discrete))
    if !(isnothing(cutoff_inf))
        bins_discrete = bins_discrete[1:cutoff_inf-1]
        counts_discrete = counts_discrete[1:cutoff_inf-1]
    end

    #cutoff_small
    cutoff_small = -5
    cutoff_small_idx = findfirst(c -> c < cutoff_small,log10.(counts_discrete))
    #log10.(counts_discrete)[cutoff_small_idx]
    if !(isnothing(cutoff_small_idx))
        bins_discrete = bins_discrete[1:cutoff_small_idx-1]
        counts_discrete = counts_discrete[1:cutoff_small_idx-1]
    end

    @info "Number of bins: $(length(bins_discrete))"

    #question? bin spacing leads to problems?
    push!(bins,bins_discrete)
    push!(bin_counts,counts_discrete)

    #save raw bins and counts 
    bins_df = [bins_discrete...,fill(NaN,df_len-length(bins_discrete))...]
    counts_df = [counts_discrete...,fill(NaN,df_len - length(counts_discrete))...]
    distr_df[:,Symbol(dir[5:end] * "_bins")] = bins_df
    distr_df[:,Symbol(dir[5:end] * "_counts")] = counts_df

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
    elseif i == 1
        par_label = par_labels[3]
        markershape =:rect
        ms = ms_high
    else 
        par_label = ""
        markershape = markershape_default
        ms = ms_default
    end

    #---------------------fit:discrete----------------

    #discrete \tau values binning
    plot!(pl_distribs, log10.(bins_discrete), log10.(counts_discrete),
        markershape=markershape, st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.5,
        lc=color,
        alpha=1,
        ms=ms,
        label=par_label)

    #try to fit to discrete binning 
    #quadratic 
    n_τ = length(bins_discrete)
    interval = 1:n_τ #1:round(Int, 0.99 * n_τ)
    logτs = log10.(bins_discrete[interval])
    logps = log10.(counts_discrete[interval])
    ps_fit,wts_fit = fit_irls(quadratic, logτs, logps, param_guesses;maxiter=maxiter,abstol=abstol)

    #=
    #prob from CurveFit
    prob = NonlinearCurveFitProblem(quadratic, param_guesses, logτs, logps)
    sol = solve(prob)
    push!(sols,sol)
    ps_fit = sol.u
    R2 = R_squared(sol.resid, logps)#calc R-squared 
    resid_ss = rss(sol;weighted=false)
    =#

    #add fit data to dataframe 
    q1,q2 = ps_fit[[2,3]]
    nr_bins = length(bins_discrete)
    push!(fitparams_df,(data_dirs[i],Δt_chart,nr_bins,min_lt,max_lt,nr_tot_lts,q1,q2,q2/q1))

    #linfit 
    interval_lin = 3:round(Int, 0.99 * n_τ)
    logτs_lin = log10.(bins_discrete[interval_lin])
    logps_lin = log10.(counts_discrete[interval_lin])
    prob_lin = CurveFitProblem(logτs_lin, logps_lin)
    sol_lin = solve(prob_lin, LinearCurveFitAlgorithm())
    ps_fit_lin = sol_lin.u
    R2_lin = R_squared(sol_lin.resid, logps_lin)#calc R-squared 
    resid_ss_lin = rss(sol_lin;weighted=false)

    if i == 2 || i == 14
        plot!(pl_distribs_lin, log10.(bins_discrete), log10.(counts_discrete),
        markershape=markershape, #st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.5,
        lc=color,
        alpha=1,
        ms=ms,
        label=par_label)
    end

    #plot fits
    τs_fit = log10.(bins_discrete[interval[1]]:0.1:bins_discrete[interval[end]])
    plot!(pl_distribs,τs_fit, quadratic(τs_fit,ps_fit),ls=:dot, lc=color, lw=2, label="")


    #plot just two fits
    if i == 1 || i == length(data_dirs)-1
        plot!(pl_distribs_part, log10.(bins_discrete), log10.(counts_discrete),
            markershape=markershape, st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            markerstrokewidth=0.5,
            lc=color,
            alpha=1,
            ms=ms,
            label=par_label)
        
        plot!(pl_distribs_part,τs_fit, quadratic(τs_fit,ps_fit),ls=:dot, lc=color, lw=2, label="")
    end

    #=
    τs_fit_lin = log10.(bins_discrete[interval_lin[1]]:0.1:bins_discrete[interval_lin[end]])
    if i == 2 || i == 14
        plot!(pl_distribs_lin,τs_fit_lin, ps_fit_lin[1] .* τs_fit_lin .+ ps_fit_lin[2],ls=:dot, lc=color, lw=2, label="")
    end

        push!(q2s,ps_fit[end])
    push!(q1s,ps_fit[2])
    push!(fit_params,ps_fit)
    push!(R2_vals,R2)
    push!(resid_ss_vals,resid_ss)
    push!(R2_lin_vals,R2_lin)
    push!(resid_ss_lin_vals,resid_ss_lin)
    @show ps_fit
    =#



    #=
    #---------------------fit:exp----------------

    #try using exponential binning
    #get exponential binning automatically adjusted for each
    bins_exp, counts_exp = autohist(lifetimes; m=1, nbins=nbins_exp,mode=:pdf)

    #exp binning
    @show length(bins_exp[1:end-1]),length(counts_exp)
    plot!(pl_distribs_exp, log10.(bins_exp[1:end-1]), log10.(counts_exp),
        markershape=markershape, #st=:scatter,
        mc=color, markerstrokecolor=:gray10,
        markerstrokewidth=0.5,
        lc=color,
        alpha=1,
        ms=ms,
        label=par_label)

    #try to fit to discrete binning 
    #quadratic 
    n_τ = length(bins_exp)
    interval = 1:round(Int, 0.99 * n_τ)
    param_guesses = [-1.0, 1.0, -1.0]
    prob_exp = NonlinearCurveFitProblem(quadratic, param_guesses, log10.(bins_exp[1:end-1]), log10.(counts_exp))
    sol_exp = solve(prob_exp)
    ps_fit_exp = sol_exp.u

    #plot fits
    xs_fit = log10.(bins_exp[1]:0.1:bins_exp[end])
    plot!(pl_distribs_exp,xs_fit, quadratic(ps_fit_exp, xs_fit),ls=:dot, lc=color, lw=2, label="")

    push!(q2s_exp,ps_fit_exp[end])
    push!(q1s_exp,ps_fit_exp[2])
    push!(fit_params_exp,ps_fit_exp)
    =#

end

pl_distribs 
pl_distribs_part 

CSV.write("../../data/distrib_data_14rval_deltat2_raw_bins_counts.csv",distr_df)

savefig(pl_distribs,figs_dir*"lifetime_distributions_and quad_params_B2_to_B4_rcurve_charting_time_$(Δt_chart).png")
savefig(pl_distribs_part,figs_dir*"lifetime_distributions_and quad_params_only_B2_and_B4_rcurve_charting_time_$(Δt_chart).png")

#plot fit goodness
pl_fitgood = plot(rs[1:(end-1)], R2_vals[1:(end-1)], markershape=:circle, 
    lc=:gray10, ylims=(0.5, 1.05), mc=:gray10, label="", xlabel=L"r", 
    ylabel=L"R^2", guidefontsize=25, left_margin=5Plots.mm, 
    xticks=[4.0:2:16;],tickfontsize=12)

pl_RSS = plot(rs[1:(end-1)], resid_ss_vals[1:(end-1)], markershape=:circle, 
    lc=:gray10, mc=:gray10, label="par", xlabel=L"r", 
    ylabel=L"RSS", guidefontsize=25, left_margin=5Plots.mm, 
    xticks=[4.0:2:16;],tickfontsize=12)

hline!(pl_fitgood, [1.0], ls=:dash, lc=:red, lw=2, label="")

#plot fit goodness
pl_fitgood_lin = plot(rs[1:(end-1)], R2_lin_vals[1:(end-1)], markershape=:circle, 
    lc=:gray10, ylims=(0.1, 1.05), mc=:gray10, label="par", xlabel=L"r", 
    ylabel=L"R^2", guidefontsize=25, left_margin=5Plots.mm, 
    xticks=[4.0:2:16;],tickfontsize=12)

pl_RSS_lin = plot(rs[1:(end-1)], resid_ss_lin_vals[1:(end-1)], markershape=:circle, 
    lc=:gray10, mc=:gray10, label="", xlabel=L"r", 
    ylabel=L"RSS", guidefontsize=25, left_margin=5Plots.mm, 
    xticks=[4.0:2:16;],tickfontsize=12)

plot!(pl_RSS,rs[1:(end-1)], resid_ss_lin_vals[1:(end-1)], markershape=:circle, 
    lc=:gray50, mc=:gray50, label="lin", xlabel=L"r", 
    ylabel=L"RSS", guidefontsize=25, left_margin=5Plots.mm, 
    xticks=[4.0:2:16;],tickfontsize=12)

hline!(pl_fitgood_lin, [1.0], ls=:dash, lc=:red, lw=2, label="")


#plot quadratic coefficients in terms of r
pl_bfit_r = plot(rs[1:end-1], abs.(q2s[1:end-1]), markershape=:circle, lc=:gray10, ms=3, label="",
    mc=:gray10, xlabel=L"r", ylabel=L"|q_2|",ylims=(-5,22)) #yaxis=:log10, xaxis=:log10)


#overlay markers of B2, B3, B4 
plot_B3 = true
plot!([rs[2]], [abs.(q2s[2])], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label=L"B_2")
plot_B3 && plot!([rs[1]], [abs.(q2s[1])], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label=L"B_3")
plot!([rs[end-4]], [abs.(q2s[end-4])], st=:scatter,markershape = :rect, mc=:purple, ms=8, label=L"B_4")


#plot quadratic coefficients in terms of w
c0 = 2.4;α0 = 0.005
wg_pairs = rescaled_params.(rs,rs,1.0,c0,α0,N)
ws = [wg_pairs[i][1] for i in 1:length(wg_pairs)]
pl_bfit_w = plot(ws[1:end-1], abs.(q2s[1:end-1]), markershape=:circle, lc=:gray10, ms=3, 
    mc=:gray10, xlabel=L"w", ylabel=L"|q_2|", ylims=(-5,50),label="") #yaxis=:log10, xaxis=:log10)

##overlay markers of B2, B3, B4 
plot_B3 = false
plot!([ws[2]], [abs.(q2s[2])], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label=L"B_2")
plot_B3 && plot!([ws[1]], [abs.(q2s[1])], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label=L"B_3")
plot!([ws[end-4]], [abs.(q2s[end-4])], st=:scatter,markershape = :rect, mc=:purple, ms=8, label=L"B_4")

#extract fitparams
q1s = fitparams_df.q1
q2s = fitparams_df.q2
q_ratios = fitparams_df.q1

ms_fitp = 5
#plot quadratic coefficients in terms of r * \Delta t (charting time)
pl_bfit_Tx = plot(rs .* Δt_chart, q2s, lc=:gray10, lw=1,
    xlabel=L"R_\Delta", ylabel=L"q_2",label="",legend=true,
    dpi=300,xlims=(6,32)) #only line

#lines for B2, B3, B4 
vline!([rs[end] * Δt_chart],lc=:orange,ls=:dash,label="",lw=2)
vline!([rs[end-1] * Δt_chart],lc=:lime,ls=:dash,label="",lw=2)
vline!([rs[1] * Δt_chart],lc=:purple,ls=:dash,label="",lw=2)

plot!(pl_bfit_Tx,rs .* Δt_chart, q2s, markershape=:circle, st=:scatter, ms=ms_fitp, 
    mc=:gray10,label=L"q_2",) #xlabel=L"Δt / T_x" 

#linear parameter 
lincolor=:gray70
plot!(pl_bfit_Tx,rs .* Δt_chart, q1s, lc=lincolor,label="")
plot!(pl_bfit_Tx,rs .* Δt_chart, q1s,markerstrokecolor=:gray10, 
    st=:scatter, markershape=:rect, mc=lincolor, ylabel=L"q_1,q_2", ms=ms_fitp,
    xlabel=L"R_\Delta", guidefontsize=25,label=L"q_1 - 1",lw=1,legendfontsize=15)




#=
#overlay markers of B2, B3, B4 
plot_B3 = true
plot!([rs[end-1]] * Δt_chart, [q2s[end-1]], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label="")
plot_B3 && plot!([rs[end]] * Δt_chart, [q2s[end]], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label="")
plot!([rs[1] * Δt_chart], [q2s[1]], st=:scatter,markershape = :rect, mc=:purple, ms=8, label="")

plot_B3 = true
plot!([rs[end-1]] * Δt_chart, [q1s[end-1]], st=:scatter,markershape=:dtriangle, mc=:lime, ms=8, label="")
plot_B3 && plot!([rs[end]] * Δt_chart, [q1s[end]], st=:scatter,markershape = :utriangle, mc=:orange, ms=8, label="")
plot!([rs[1] * Δt_chart], [q1s[1]], st=:scatter,markershape = :rect, mc=:purple, ms=8, label="")
=#

savefig(pl_bfit_Tx,figs_dir*"fit_parameters_q1_q2_B2_to_B4_rcurve_charting_time_$(Δt_chart).png")



#=
#add inset to left panel 
plot!(pl_bfit_Tx,rs .* Δt_chart, abs.(q2s) ./ q1s, 
    inset = (1, bbox(0.5,0.05,0.5,0.5)), subplot = 2, #ylims=(-3,3),
    lc=:gray10, markershape=:circle, mc=:gray10, ylabel=L"|q_2|/q_1", 
    xlabel=L"R_\Delta", guidefontsize=25,label="",lw=1)
=#

plot(rs .* Δt_chart, q_ratios, 
    lc=:gray10, markershape=:circle, mc=:gray10, ylabel=L"q_2/q_1", 
    xlabel=L"R_\Delta", guidefontsize=25,label="",lw=1,ylims=(-5,5))

plot(rs .* Δt_chart, q1s, 
    lc=:gray10, markershape=:circle, mc=:gray10, ylabel=L"q_1", 
    xlabel=L"R_\Delta", guidefontsize=25,label="",lw=1,ylims=(-5,5))

#ylims=(2,9), yticks=[3,6,9],xticks=[0.5,0.35,0.2],
#xlims=(0.15,0.5)

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
annotate!(pl_distfit_sup[1], 1.71, 0.5, text(L"r", :gray10, :center, annot_size)) #colorbar label

savefig(pl_distfit_sup,figs_dir*"lifetime_distributions_and quad_params_B2_to_B4_change_r_charting_time_$(Δt_chart).png")