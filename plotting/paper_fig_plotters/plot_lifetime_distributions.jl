using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#read lifetime filenames
files = readdir(data_dir_lt)
files_lifetime = files[findall(f -> occursin(".dat", f), files)]#[[13,14,15,16,1:12...]] #select binary only
files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times] #select binary only

ms = 5
linealphas = [0.3,0.3,0.3,1.0] 
plots_loglog = [] #gather plots here
plots_linhist = [] #gather plots here
plots_linlin = [] #gather plots here
plots_loglin = [] #gather plots here
lifetimes_stats = []
#plot_order = [1,2,3,4,5]
#plot_fits = false
#fileset_order = [1:4;]

#fit curves
#quadratic(a, x) = @. a[1] + a[2] * x  + a[3] * x^2 #CurveFit
quadratic(x,p) = @. p[1] + p[2] * x  + p[3] * x^2 #LsqFit
cubic(a, x) = @. a[1] + a[2] * x^3
fitparam_dataframes = []
param_guesses_quad = [-10.0, 1.0, -10.0]
param_guesses_lin = [1.0,-1.0]

#loop through lifetime files sorted by charting time
for (i, fileset) in enumerate(files_lifetime_by_time)

    pls_loglog = [plot(), plot(), plot()]#init empty 3-plot
    pls_linlin = [plot(), plot(), plot()]#init empty 3-plot
    pls_loglin = [plot(), plot(), plot()]#init empty 3-plot

    #dataframes 
    df = DataFrame(parameter=String[],chart_time=Float64[],q1=Float64[],q2=Float64[],q_ratio=Float64[])

    ylabel = i == 1 ? L"\log(p(\tau))" : ""

    chart_time = chart_times[i]
    println("Plots for charting time: $(chart_time)")

    #take every file with a given charting time
    for (j, f) in enumerate(fileset[fileset_order])
        println("File: " * f[1:22] * "...")
        #read lifetimes
        lifetimes = reinterpret(Int64, read(data_dir_lt * f)) 
        color = colors[j]

        #--------------------------get binning for lifetimes--------------------------

        #preset bins, one for each possible \tau value (discrete)
        min_lt,max_lt = extrema(lifetimes)
        bin_edges = min_lt:1:max_lt

        #plot normal linear histogram plots
        pl_hist = barhist(lifetimes; bins=100,
            label="$(par_labels[j])",
            title=L"\Delta t_L = %$(chart_times[i])",
            color=color)
        push!(plots_linhist, pl_hist)

        #skip if distribution collapses into one or two bins 
        length(bin_edges) < 2 && continue

        #lifetime stats
        nr_lifetimes =  length(lifetimes)

        push!(lifetimes_stats, ["$(par_labels[j])", chart_time,max_lt, nr_lifetimes])

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

        
        #---------------exponential binning----------------

        #get exponential binning automatically adjusted for each
        bins_exp, counts_exp = autohist(lifetimes; m=1, nbins=nbins,mode=:pdf)
    
        #---------------linear binning (coarse)----------------

        #get linear binning automatically adjusted for each
        h_lin = fit(Histogram, lifetimes; nbins=nbins)
        h_lin = normalize(h_lin;mode=:pdf)
        bins_lin = h_lin.edges[1][1:end-1]
        counts_lin = h_lin.weights
        

        #--------------------------------make plots-------------------------------------------

        #--------------------plot uniform discrete binning-------------------------- 

        plot!(pls_loglog[1], log10.(bins_discrete), log10.(counts_discrete), 
            markershape=markershapes[j],st=:scatter,
            mc=color, 
            markerstrokecolor=:gray10,
            markerstrokewidth=msw,
            lc=color,
            alpha = 1,
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t_L = %$(chart_time)", 
            ms=ms,
            label= par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        n_τ = length(bins_discrete)

        #add fits, guiding lines to the largest chart times (i==3)

        #add linear fit, with offset (on the side)
        #only for the largest charting time (i == 3) 
        if plot_fits && i ∈ chart_time_fit_idx && j == linfit_idx

            #linear 
            linear_interval = round(Int, 0.2 * n_τ):round(Int, 0.8 * n_τ)

            b, a = linreg(log10.(bins_discrete[linear_interval]) .+ 0.2, log10.(counts_discrete[linear_interval]))
            @info "Parameters from linear fit: $(b),$(a)" 
            a_round = round(a;digits=3)

            plot!(log10.(bins_discrete[linear_interval]), log10.(10 .^ (b) .* bins_discrete[linear_interval] .^ (-4)),
                ls=:dash, lc=:gray10, lw=1, label="", legendfontsize=15)

        end

        #quadratic
        if plot_fits && i ∈ chart_time_fit_idx && j ∈ quadfit_idxs
            
            #quadratic 
            interval = 1:round(Int,  n_τ)
            
            #=
            #nonadaptive simple least squares fit
            prob = NonlinearCurveFitProblem(quadratic, param_guesses, log10.(bins_discrete[interval]), log10.(counts_discrete[interval]))
            sol = solve(prob)
            ps_fit = sol.u
            =#

            logτs = log10.(bins_discrete[interval])
            logps = log10.(counts_discrete[interval])
            ps_fit,wts_fit = fit_irls(quadratic, logτs, logps, param_guesses_quad;maxiter=maxiter,abstol=abstol)

            @info "Parameters from parabolic fit (q1,q2): $(ps_fit[2]),$(ps_fit[3])"
            logτs_fit = logτs[1]:0.01:logτs[end]
            plot!(logτs_fit, quadratic(logτs_fit,ps_fit),
                ls=:dot, lc=color, lw=2, label="", legendfontsize=15)

            #save fitparams 
            push!(df,(par_labels_string[j],chart_time,ps_fit[2],ps_fit[3],ps_fit[3]/ps_fit[2]))

        end

        #=
        #add higher order fit
        if plot_fits && i == 1 && j == 2
            interval = 1:round(Int, 0.02 * n_τ)
            #vline!([log10.(bins_discrete[quad_interval[1]]), log10.(bins_discrete[quad_interval[end]])], ls=:dash)
            prob = NonlinearCurveFitProblem(cubic, param_guesses, log10.(bins_discrete[interval]), log10.(counts_discrete[interval]))
            sol = solve(prob)
            ps_fit = sol.u
            plot!(log10.(bins_discrete[interval]), cubic(ps_fit,log10.(bins_discrete[interval])),
                ls=:dot, lc=:gray10, lw=2, label="", legendfontsize=15)
        end
        =#

        plot!(pls_loglin[1], bins_discrete, log10.(counts_discrete),
            markershape=markershapes[j],st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_linlin[1], bins_discrete, counts_discrete,
            markershape=markershapes[j],st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)
        
        #show bin sizes 
        #plot!(pls[1],log10.(bins_discrete), fill(-11-j,length(bins_discrete)),
        #    markershape=:vline, st=:scatter,mc=color,ms=5,
        #    markerstrokewidth = 0.5,label="")


        #--------------------plot exponential binning-------------------------- 

        plot!(pls_loglog[2], log10.(bins_exp[1:end-1]), log10.(counts_exp),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        plot!(pls_linlin[2], bins_exp[1:end-1], counts_exp,
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_loglin[2], bins_exp[1:end-1], log10.(counts_exp),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        #show bin sizes 
        #plot!(pls[2], log10.(bins_exp), fill(-11-j, length(bins_exp)),
        #    markershape=:vline, st=:scatter, mc=color, ms=5,
        #    markerstrokewidth=0.5, label="")

        #---------------------plot linear binning (coarse)--------------------

        plot!(pls_loglog[3], log10.(bins_lin), log10.(counts_lin), 
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        plot!(pls_linlin[3], bins_lin, counts_lin,
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_loglin[3], bins_lin, log10.(counts_lin),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            lc=color,
            dpi=300,
            title=L"\Delta t_L = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        #show bin sizes 
        #plot!(pls[3], log10.(bins_lin), fill(-11-j, length(bins_lin)),
        #    markershape=:vline, st=:scatter, mc=color, ms=5,
        #    markerstrokewidth=0.5, label="")

        
    end
    push!(plots_loglog,pls_loglog)
    push!(plots_linlin, pls_linlin)
    push!(plots_loglin, pls_loglin)
    push!(fitparam_dataframes,df)
end

n_ct = length(chart_times)

#---------------------------loglog-------------------------
plots_reorder_loglog = [plots_loglog[i][j] for j in 1:3 for i in 1:n_ct]
xlimits = (0.0, 3.2); ylimits = (-12.0, 0.0);
pl_loglog_all = plot(plots_reorder_loglog..., layout=(3, n_ct), size=(1200, 800),xlims=xlimits,#ylims=ylimits,
    plot_title=L"N=%$N,k=%$k", left_margin=7Plots.mm,legend=:topright,legendfontsize=8,dpi=300)

#---------------------------linlin-------------------------
plots_reorder_linlin = [plots_linlin[i][j] for j in 1:3 for i in 1:n_ct]
xlimits = (0.0, 3.2);ylimits = (-12.0, 0.0);
pl_linlin_all = plot(plots_reorder_linlin..., layout=(3, n_ct), size=(1200, 800),xlims=(0,180),
    plot_title=L"N=%$N,k=%$k", left_margin=5Plots.mm, legend=:topright, legendfontsize=8, dpi=300)

#---------------------------loglin-------------------------
plots_reorder_loglin = [plots_loglin[i][j] for j in 1:3 for i in 1:n_ct]
xlimits = (0.0, 3.2);
ylimits = (-12.0, 0.0);
pl_loglin_all = plot(plots_reorder_loglin..., layout=(3, n_ct), size=(1200, 800), xlims=(0, 800),
    plot_title=L"N=%$N,k=%$k", left_margin=5Plots.mm, legend=:topright, legendfontsize=8, dpi=300)


#savedir_pl = "../../paper_figs/lifetime_distributions_rescaled_horizontal_vertical/"
#mkpath(savedir_pl)
#param_string = files[2][20:end-16]
param_string = files[2][56:end-16]

#savefig(pl_loglog_all, savedir_pl * "lifetime_distributions_integrated_" * param_string * "_log_log_discrete_exp_linbin.png")
#savefig(pl_linlin_all, savedir_pl * "lifetime_distributions_integrated_" * param_string * "_lin_lin_discrete_exp_linbin.png")
#savefig(pl_loglin_all, savedir_pl * "lifetime_distributions_integrated_" * param_string * "_log_lin_discrete_exp_linbin.png")

