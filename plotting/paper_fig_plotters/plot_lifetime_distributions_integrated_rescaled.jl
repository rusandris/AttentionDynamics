using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("../../scripts/new/chart_functions.jl")

#general parameters
#colors = [:lime, :cyan, :orange, :darkgreen]
colors = [:lime, :cyan, :gray10]
markershapes = [:circle, :cross]
nbins = 30 #[25,10]

#lifetimes path 
data_dir_traj_root = "../../data/lifetimes_precomp_rescaled_vertical/"
#read params
files = readdir(data_dir_traj_root)
files_lifetime = files[findall(f -> occursin(".dat", f), files)] #select binary only
params = readdlm(data_dir_traj_root * files[1])
#params = [0.45 22.5; 0.5 22.5; 0.55 22.5]
ws = params[:,1]
gs = params[:,2]
gs_rounded = round.(gs; digits=4)
ws_rounded = round.(ws; digits=4)


#traj specific parameters (not from filename, hardcoded)
N = 10
Δt_traj = 0.1 #sampling of trajectories
T = 100000
t = [0.0:Δt_traj:T;]
data_shape = (Int(T / Δt_traj) + 1, N)
k = 4

plots = []
#pls = [plot(), plot(), plot()]

chart_times = [0.1, 0.5 ,1.0] #[1.0,2.0]
#chart_time = 0.1 #0.1 #1.0 #0.1 #0.5
files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times]#select binary only


for (i, fileset) in enumerate(files_lifetime_by_time)
    pls = [plot(), plot(), plot()]#init empty 3-plot
    chart_time = chart_times[i]
    for (j,f) in enumerate(fileset)

        #read lifetimes
        lifetimes = reinterpret(Int64, read(data_dir_traj_root * f))
        color = colors[j]

        #--------------------------get binning for lifetimes--------------------------

        #preset bins, one for each possible \tau value (discrete)
        bin_edges = minimum(lifetimes):1:maximum(lifetimes)
        h_discrete = fit(Histogram, Float64.(lifetimes), bin_edges)
        bins_discrete = h_discrete.edges[1][1:end-1]
        counts_discrete = h_discrete.weights

        #get exponential binning automatically adjusted for each
        bins_exp, counts_exp = autohist(lifetimes; m=1, nbins=nbins)
        #normalize by bin width
        counts_exp = counts_exp ./ diff(bins_exp)
        bins_exp = sqrt.([bins_exp[i] * bins_exp[i+1] for i in 1:length(bins_exp)-1]) #right edges

        #get linear binning automatically adjusted for each
        h_lin = fit(Histogram, lifetimes; nbins=nbins)

        bins_lin = h_lin.edges[1]
        counts_lin = h_lin.weights
        
        counts_lin = counts_lin ./ diff(bins_lin)
        bins_lin = sqrt.([bins_lin[i] * bins_lin[i+1] for i in 1:length(bins_lin)-1]) #right edges

        #make plots 
        plot!(pls[1], log10.(bins_discrete), log10.(counts_discrete), 
            markershape=:circle,#st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=4,
            label=L"w = %$(ws_rounded[i]), g = %$(gs_rounded[i])",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot!(pls[2], log10.(bins_exp), log10.(counts_exp),
            markershape=:circle, #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=4,
            label=L"w = %$(ws_rounded[i]), g = %$(gs_rounded[i])",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot!(pls[3], log10.(bins_lin), log10.(counts_lin), 
            markershape=:circle, #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=4,
            label=L"w = %$(ws_rounded[i]), g = %$(gs_rounded[i])",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")
    end
        push!(plots,pls)
end


plots_reorder = [plots[i][j] for j in 1:3 for i in 1:3]
pl_all = plot(plots_reorder..., layout=(3, 3), size=(1000, 800),
    plot_title=L"N=%$N,k=%$k", dpi=300)

savedir = "../../paper_figs/lifetime_distributions_rescaled/"
mkpath(savedir)
pl_all = plot(pls...,layout=(3,1),size=(400,800))
savefig(pl_all, savedir * "lifetime_distributions_rescaled_k$(k)_change_orig_param_charting time_$(chart_time)_integrated_log_log_discrete_exp_linbin_verticalB_10ic.png")


#=
bin_edges = minimum(lifetimes):1:maximum(lifetimes)
h = fit(Histogram,Float64.(lifetimes),bin_edges)

hist_lin = plot(h.edges[1][1:end-1],h.weights,mc=:gray10,ms=6,xlims=(0,25),
           ylabel="count",xlabel=L"\tau",guidefontsize=18,st=:scatter,
           label=L"w=%$(ws_rounded[1]),g=%$(gs_rounded[1])",legendfontsize=12,legend=:topleft)

plot!(h2.edges[1][1:end-1],h2.weights,mc=:orange,ms=6,xlims=(0,25),
           ylabel="count",xlabel=L"\tau",guidefontsize=18,st=:scatter,
           label=L"w=%$(ws_rounded[2]),g=%$(gs_rounded[2])",
           markerstrokecolor=:orange,ylims=(0,8*1e5))

hist_log = plot(log10.(h.edges[1][1:end-1]),log10.(h.weights),mc=:gray10,ms=6,
        ylabel="log(count)",xlabel=L"\log(\tau)",guidefontsize=18,st=:scatter,
        label=L"w=%$(ws_rounded[1]),g=%$(gs_rounded[1])",legendfontsize=12,legend=:topright)

plot!(log10.(h2.edges[1][1:end-1]),log10.(h2.weights),mc=:orange,ms=6,
    ylabel="log(count)",xlabel=L"\log(\tau)",guidefontsize=18,st=:scatter,
    label=L"w=%$(ws_rounded[2]),g=%$(gs_rounded[2])",legendfontsize=12,
    legend=:topright,markerstrokecolor=:orange)

pl = plot(hist_lin, hist_log, layout=(1, 2), size=(1200, 400), dpi=300, bottom_margin=6Plots.mm, left_margin=6Plots.mm)=#

#=

colors = [:purple, :red, :orange]

for (i,f) in enumerate(files_integ[[6,10,2]])
    chart_time = chart_times[i]
    lifetimes_c_integ = vec(readdlm(data_dir_lt * "/integrated/" * f))
    color = colors[i]
    plot_distrib!(pl, lifetimes_c_integ, chart_time * Δt;
        color=color, nbins=nbins, mode=:none,
        scale=:log, binning=:exp,st=:path,
        title=L"k=20", ms=4,
        label_par=L"c = %$(cs[2]), Δt=%$(chart_time*Δt)",
        xlabel=L"log(\tau \cdot \Delta t)")
end
plot!(pl,legend=:outertopright)
=#

