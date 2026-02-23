using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
cd(@__DIR__)

#general parameters
#colors = [:lime, :cyan, :orange, :darkgreen]
#colors = [:lime, :cyan, :orange, :gray40]
#colors = [:lime, :gray10, :cyan, :purple, :blue]
colors = [:lime, :cyan, :orange]
markershapes = [:circle, :cross]
nbins = 50 #[25,10]

#lifetimes path
#data_dir_lt_root = "../../data/lifetimes_precomp_orig/"
#data_dir_lt_root = "../../data/lifetimes_precomp_rescaled/"

files_traj = readdir(data_dir_traj)


#p_values = [params[:, i] for i in 1:size(params)[2]]
p_values = [[params_changed...]] #[gs[1:3],ws[1:3]]
p_names = ["c"] #["g","w"] #["rp", "rc"]
p_names_ltx = [L"c"] #["g", "w"] #[L"r_p", L"r_c"]


#traj specific parameters 
N, Δt_traj,T = readdlm(data_dir_traj*"traj_spec_params.txt")
N = Int(N)
t = [0.0:Δt_traj:T;]
data_shape = (Int(T / Δt_traj) + 1, N)

plots = []
#read chart times 
chart_times = readdlm(data_dir_lt * "chart_times.txt")

#read lifetime filenames
files = readdir(data_dir_lt)
files_lifetime = files[findall(f -> occursin(".dat", f), files)]#[[13,14,15,16,1:12...]] #select binary only
files_lifetime_by_time = [files_lifetime[findall(f -> occursin("charting_$time", f), files_lifetime)] for time in chart_times]#select binary only

#label names  

par_labels = []
for j in 1:length(p_values[1])
    #put strings together
    par_names = ""
    for i in 1:length(p_values)
        par_names *= p_names_ltx[i] * L" = %$(p_values[i][j])"
    end
    @show par_names
    push!(par_labels,par_names)
end

for (i, fileset) in enumerate(files_lifetime_by_time[1:3])
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
            label= par_labels[j],
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        #show bin sizes 
        plot!(pls[1],log10.(bins_discrete), fill(-0.5*j,length(bins_discrete)),
            markershape=:vline, st=:scatter,mc=color,ms=5,
            markerstrokewidth = 0.5,label="")

        plot!(pls[2], log10.(bins_exp), log10.(counts_exp),
            markershape=:circle, #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=4,
            label=p_names_ltx[1] * L" = %$(p_values[1][j])," ,#* p_names_ltx[2] * L" = %$(p_values[2][j])",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count/\Delta b_i)")

        #show bin sizes 
        plot!(pls[2], log10.(bins_exp), fill(-0.5*j, length(bins_exp)),
            markershape=:vline, st=:scatter, mc=color, ms=5,
            markerstrokewidth=0.5, label="")

        plot!(pls[3], log10.(bins_lin), log10.(counts_lin), 
            markershape=:circle, #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"\Delta t = %$(chart_time)", ms=4,
            label=p_names_ltx[1] * L" = %$(p_values[1][j]),", # * p_names_ltx[2] * L" = %$(p_values[2][j])",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count/\Delta b_i)")

        #show bin sizes 
        plot!(pls[3], log10.(bins_lin), fill(-0.5*j, length(bins_lin)),
            markershape=:vline, st=:scatter, mc=color, ms=5,
            markerstrokewidth=0.5, label="")
    end
        push!(plots,pls)
end

plots_reorder = [plots[i][j] for j in 1:3 for i in 1:3]
xlimits = (0.0, 3.2); ylimits = (-2.0, 7.0);
pl_all = plot(plots_reorder..., layout=(3, 3), size=(1200, 800),xlims=xlimits,ylims=ylimits,
    plot_title=L"N=%$N,k=%$k", left_margin=5Plots.mm,legend=:topright,legendfontsize=5,dpi=300)

#savedir_pl = "../../paper_figs/lifetime_distributions_rescaled_horizontal_vertical/"
#mkpath(savedir_pl)
#param_string = files[2][20:end-16]
param_string = files[2][25:end-16]
savedir_pl = savedir_pl_root * "c_k$(k)_Δt$(Δt_traj)_T$(T)/"
mkpath(savedir_pl)
savefig(pl_all, savedir_pl * "lifetime_distributions_integrated_" * param_string * "_log_log_discrete_exp_linbin_change_c_ic$(nr_ics).png")


