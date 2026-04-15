using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#read lifetime filenames
files = readdir(data_dir_lt)
files_lifetime = files[findall(f -> occursin(".dat", f), files)]#[[13,14,15,16,1:12...]] #select binary only
files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times]#select binary only

ms = 5
αs = [0.3,0.3,0.3,0.3,1.0] 
plots_loglog = [] #gather plots here
plots_linhist = [] #gather plots here
plots_linlin = [] #gather plots here
plots_loglin = [] #gather plots here

for (i, fileset) in enumerate(files_lifetime_by_time)

    pls_loglog = [plot(), plot(), plot()]#init empty 3-plot
    pls_linlin = [plot(), plot(), plot()]#init empty 3-plot
    pls_loglin = [plot(), plot(), plot()]#init empty 3-plot

    ylabel = i == 1 ? L"\log(\textnormal{PDF})" : ""

    chart_time = chart_times[i]
    @show chart_time
    for (j,f) in enumerate(fileset[[1,3,4,5,2]])
        @show f
        #read lifetimes
        lifetimes = reinterpret(Int64, read(data_dir_lt * f)) 
        @show length(lifetimes)
        color = colors[j]

        #--------------------------get binning for lifetimes--------------------------

        #preset bins, one for each possible \tau value (discrete)
        min_lt,max_lt = extrema(lifetimes)
        bin_edges = min_lt:1:max_lt

        #plot normal linear histogram plots
        pl_hist = barhist(lifetimes; bins=100,
            label="$(par_labels[j])",
            title=L"\Delta t = %$(chart_times[i])",
            color=color)
        push!(plots_linhist, pl_hist)

        #skip if distribution collapses into one or two bins 
        length(bin_edges) <= 2 && continue

        h_discrete = fit(Histogram, Float64.(lifetimes), bin_edges)
        @show extrema(lifetimes)
        h_discrete = normalize(h_discrete; mode=:pdf)
        bins_discrete = h_discrete.edges[1][1:end-1]
        counts_discrete = h_discrete.weights

        
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
            markershape=markershapes[j],#st=:scatter,
            mc=color, markerstrokecolor=:gray10,
            markerstrokewidth=0.1,
            lc=color,
            alpha = αs[j],
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t = %$(chart_time)", 
            ms=ms,
            label= par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        #add linear guiding line
        if i == 1 && j == 5
            n_τ = length(bins_discrete)
            linear_interval = round(Int, 0.05 * n_τ):round(Int, 0.15 * n_τ)

            #vline!([log10.(bins_discrete[linear_interval[1]]), log10.(bins_discrete[linear_interval[end]])],ls=:dash)
            b, a = linreg(log10.(bins_discrete[linear_interval]), log10.(counts_discrete[linear_interval]))
            @show b,a

            plot!(log10.(bins_discrete[linear_interval]), log10.(10 .^ (b+0.3) .* bins_discrete[linear_interval] .^ (-4)),
                ls=:solid, lc=:red,lw=1.5, label=L"\tau ^{-4}", legendfontsize=15)
        end

        plot!(pls_loglin[1], bins_discrete, log10.(counts_discrete),
            markershape=markershapes[j],#st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_linlin[1], bins_discrete, counts_discrete,
            markershape=markershapes[j],#st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"N = %$(N), k = %$(k), \Delta t = %$(chart_time)", ms=ms,
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
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        plot!(pls_linlin[2], bins_exp[1:end-1], counts_exp,
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_loglin[2], bins_exp[1:end-1], log10.(counts_exp),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
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
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=ylabel)

        plot!(pls_linlin[3], bins_lin, counts_lin,
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\tau",
            ylabel=ylabel)

        plot!(pls_loglin[3], bins_lin, log10.(counts_lin),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
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

end

n_ct = length(chart_times)

#---------------------------loglog-------------------------
plots_reorder_loglog = [plots_loglog[i][j] for j in 1:3 for i in 1:n_ct]
xlimits = (0.0, 3.2); ylimits = (-12.0, 0.0);
pl_loglog_all = plot(plots_reorder_loglog..., layout=(3, n_ct), size=(1200, 800),xlims=xlimits,#ylims=ylimits,
    plot_title=L"N=%$N,k=%$k", left_margin=5Plots.mm,legend=:topright,legendfontsize=8,dpi=300)

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

