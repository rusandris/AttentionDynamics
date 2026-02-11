using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("../../scripts/new/chart_functions.jl")

#general parameters
colors = [:lime, :cyan, :orange, :darkgreen]
markershapes = [:circle, :cross]

nbins = 30 #[25,10]
#chart_times = [0.1, 1.0] #[1.0,2.0]
chart_time = 0.5 #1.0 #0.1
k = 4

data_dir_traj_root = "../../data/trajs_precomp_rescaled/"
subdirs = readdir(data_dir_traj_root)
plots = []

for (dir_idx,data_dir_traj) in enumerate(subdirs)
    @show data_dir_traj
    path = data_dir_traj_root * data_dir_traj * "/"
    pls = [plot(), plot(), plot()]
    color = colors[dir_idx]

    #---------------------------read precomputed trajectories------------------------------------ 
    dir_content = readdir(path)
    files = dir_content[findall(f -> occursin(".txt",f) || occursin(".dat",f),dir_content)]

    #read params
    gs = vec(readdlm(path * files[1]))
    ws = vec(readdlm(path * files[end]))
    gs_rounded = round.(gs; digits=4)
    ws_rounded = round.(ws; digits=4)

    files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only
    nr_ics = Int(length(files_trajs)/length(gs))

    files_trajs_ics = [files_trajs[findall(f -> occursin("g$(gs[i])_w$(ws[i])", f), files_trajs)] for i in 1:length(gs)]

    #traj specific parameters from filename
    N = parse(Int64, files_trajs[1][end-34:end-33])
    Δt = parse(Float64, files_trajs[1][end-10:end-8]) #sampling of trajectories
    T = parse(Float64, files_trajs[1][end-30:end-25])
    t = [0.0:Δt:T;]
    data_shape = (Int(T / Δt) + 1, N)

    #-----------------------------compute and gather lifetimes------------------------------

    #loop through parameters (corresponding filesets)
    for (i,fileset) in enumerate(files_trajs_ics)

        lifetimes_all = Int64[]

        #read each initial condition file 
        #compute lifetimes for each
        for f in fileset
            @show f
            fname, = splitext(f)
            
            #read trajectory
            traj = Matrix(reshape(reinterpret(Float64, read(path * f)), data_shape))

            #save traj pl
            mkpath(path*"figs/")
            savefig(plot(traj[1:500,:],legend=false,title=fname,size=(1200,300)),path*"figs/"*fname*".png")
            
            #integrate activity segment to segment
            int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt))
            
            lts = get_lifetimes_integ(int_act, k)
            @show length(lts)
            push!(lifetimes_all,lts...)
        end

        #--------------------------get binning for lifetimes--------------------------

        #preset bins, one for each possible \tau value (discrete)
        bin_edges = minimum(lifetimes_all):1:maximum(lifetimes_all)
        h_discrete = fit(Histogram, Float64.(lifetimes_all), bin_edges)
        bins_discrete = h_discrete.edges[1][1:end-1]
        counts_discrete = h_discrete.weights

        #get exponential binning automatically adjusted for each
        bins_exp, counts_exp = autohist(lifetimes_all; m=1, nbins=nbins)
        #normalize by bin width
        counts_exp = counts_exp ./ diff(bins_exp)
        bins_exp = sqrt.([bins_exp[i] * bins_exp[i+1] for i in 1:length(bins_exp)-1]) #right edges

        #get linear binning automatically adjusted for each
        h_lin = fit(Histogram, lifetimes_all; nbins=nbins)

        bins_lin = h_lin.edges[1]
        counts_lin = h_lin.weights

        counts_lin = counts_lin ./ diff(bins_lin)
        bins_lin = sqrt.([bins_lin[i] * bins_lin[i+1] for i in 1:length(bins_lin)-1]) #right edges

        #make plots 
        plot!(pls[1], log10.(bins_discrete), log10.(counts_discrete), 
            markershape=markershapes[i],#st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"k=%$k", ms=4,
            label=L"w = %$(ws_rounded[i]), g = %$(gs_rounded[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot!(pls[2], log10.(bins_exp), log10.(counts_exp),
            markershape=markershapes[i], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"k=%$k", ms=4,
            label=L"w = %$(ws_rounded[i]), g = %$(gs_rounded[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot!(pls[3], log10.(bins_lin), log10.(counts_lin), 
            markershape=markershapes[i], #st=:scatter,
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"k=%$k", ms=4,
            label=L"Lw = %$(ws_rounded[i]), g = %$(gs_rounded[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")
    end
    push!(plots,pls)
end

pl_expbin = plot(plots_expbin..., layout=(1, 2), size=(800, 300), plot_title="expbin",
    legend=:bottomleft, left_margin=6Plots.mm, bottom_margin=6Plots.mm)

pl_linbin = plot(plots_linbin..., layout=(1, 2), size=(800, 300), plot_title="linbin",
    legend=:bottomleft, left_margin=6Plots.mm)

plots_reorder = [plots[i][j] for j in 1:3 for i in 1:4]
pl_all = plot(plots_reorder..., layout=(3, 4), size=(1200, 800),dpi=300)

savedir = "../../paper_figs/lifetime_distributions_rescaled/"
mkpath(savedir)
savefig(pl_all, savedir * "lifetime_distributions_rescaled_k$(k)_change_orig_param_charting time_0.5_integrated_log_log_discrete_bin_exp_bin_linbin_properbins_bin_norm_10ic.png")
savefig(pl_expbin, savedir * "lifetime_distributions_rescaled_k$(k)_change_alpha_charting time_0.1_1_integrated_log_log_exp_bin_separated_properbins_bin_norm_10ic.png")

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

