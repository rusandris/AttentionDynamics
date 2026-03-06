using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
cd(@__DIR__)

#read lifetime filenames
files = readdir(data_dir_lt)
files_lifetime = files[findall(f -> occursin(".dat", f), files)]#[[13,14,15,16,1:12...]] #select binary only
files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times]#select binary only

ms = 2
plots = [] #gather plots here
for (i, fileset) in enumerate(files_lifetime_by_time)
    pls = [plot(), plot(), plot()]#init empty 3-plot
    chart_time = chart_times[i]
    @show chart_time
    for (j,f) in enumerate(fileset)
        @show f
        #read lifetimes
        lifetimes = reinterpret(Int64, read(data_dir_lt * f)) 
        color = colors[j]

        #--------------------------get binning for lifetimes--------------------------

        #preset bins, one for each possible \tau value (discrete)
        bin_edges = minimum(lifetimes):1:maximum(lifetimes)
        h_discrete = fit(Histogram, Float64.(lifetimes), bin_edges)
        h_discrete = normalize(h_discrete; mode=:pdf)
        bins_discrete = h_discrete.edges[1][1:end-1]
        counts_discrete = h_discrete.weights

        #get exponential binning automatically adjusted for each
        bins_exp, counts_exp = autohist(lifetimes; m=1, nbins=nbins,mode=:pdf)

        #normalize by bin width
        #counts_exp = counts_exp ./ diff(bins_exp)
        #bins_exp = sqrt.([bins_exp[i] * bins_exp[i+1] for i in 1:length(bins_exp)-1]) #right edges

        #get linear binning automatically adjusted for each
        h_lin = fit(Histogram, lifetimes; nbins=nbins)
        h_lin = normalize(h_lin;mode=:pdf)
        bins_lin = h_lin.edges[1][1:end-1]
        counts_lin = h_lin.weights
        
        #counts_lin = counts_lin ./ diff(bins_lin)
        #bins_lin = sqrt.([bins_lin[i] * bins_lin[i+1] for i in 1:length(bins_lin)-1]) #right edges

        #check normalisation 
        #@show sum(diff(h_discrete.edges[1]) .* counts_discrete)
        #@show sum(diff(bins_exp) .* counts_exp)
        #@show sum(diff(h_lin.edges[1]) .* counts_lin)

        #make plots 
        plot!(pls[1], log10.(bins_discrete), log10.(counts_discrete), 
            markershape=markershapes[j],#st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label= par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=L"\log(\textnormal{pdf})")

        #show bin sizes 
        #plot!(pls[1],log10.(bins_discrete), fill(-11-j,length(bins_discrete)),
        #    markershape=:vline, st=:scatter,mc=color,ms=5,
        #    markerstrokewidth = 0.5,label="")

        plot!(pls[2], log10.(bins_exp[1:end-1]), log10.(counts_exp),
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=L"\log(\textnormal{pdf})")

        #show bin sizes 
        #plot!(pls[2], log10.(bins_exp), fill(-11-j, length(bins_exp)),
        #    markershape=:vline, st=:scatter, mc=color, ms=5,
        #    markerstrokewidth=0.5, label="")

        plot!(pls[3], log10.(bins_lin), log10.(counts_lin), 
            markershape=markershapes[j], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=ms,
            label=par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=L"\log(\textnormal{pdf})")

        #show bin sizes 
        #plot!(pls[3], log10.(bins_lin), fill(-11-j, length(bins_lin)),
        #    markershape=:vline, st=:scatter, mc=color, ms=5,
        #    markerstrokewidth=0.5, label="")
    end
    push!(plots,pls)
end

n_ct = length(chart_times)
plots_reorder = [plots[i][j] for j in 1:3 for i in 1:n_ct]
xlimits = (0.0, 3.2); ylimits = (-12.0, 0.0);
pl_all = plot(plots_reorder..., layout=(3, n_ct), size=(1200, 800),xlims=xlimits,#ylims=ylimits,
    plot_title=L"N=%$N,k=%$k", left_margin=5Plots.mm,legend=:topright,legendfontsize=8,dpi=300)

#savedir_pl = "../../paper_figs/lifetime_distributions_rescaled_horizontal_vertical/"
#mkpath(savedir_pl)
#param_string = files[2][20:end-16]
param_string = files[2][56:end-16]

savefig(pl_all, savedir_pl * "lifetime_distributions_integrated_" * param_string * "_log_log_discrete_exp_linbin.png")


