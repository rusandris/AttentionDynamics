using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")
#-------------------------plot settings------------------------------

#get lifetimes with different k, charting times 
ks = [10, 50, 100]
Δts_chart = [1, 10, 30]
savedir = "../../paper_figs/lifetime_distributions/test6/"
nbins = 50
figsize = (400, 800)
legendfontsize = 8
Δt = 1e-1 #sampling rate of trajs 

#------------------------------rp=rc------------------------------
data_dir_traj_root = "../../data/trajs_precomp_test6/"
data_dir_lt_root = "../../data/lifetimes_precomp_test6/"
data_dir_lt = data_dir_lt_root * "rp=rc/"
data_dir_traj = data_dir_traj_root * "rp=rc/"
files_trajs = readdir(data_dir_traj)

rs = vec(readdlm(data_dir_traj * files_trajs[1]))

plots_topk = []
for (k_idx, k) in enumerate(ks)

    data_files = filter(f -> occursin("Δt1.", f) && occursin("k_$(k)_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt*f)
        push!(lifetimes_x,lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = r_p = ")

    push!(plots_topk, pl_distrib_x)
end

pl_topk = plot(plots_topk..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)

savefig(pl_topk, savedir * "lifetime_distributions_r_topk_shapes_recolor_lowr.png")

#keep k the same, change charting time
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    data_files = filter(f -> occursin("Δt$(Δt_chart).txt", f) && occursin("k_10_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = r_p = ")

    push!(plots_charting, pl_distrib_x)
end

pl_ch = plot(plots_charting..., layout=(3, 1), size=figsize, legend=:bottomleft, legendfontsize=legendfontsize, dpi=300)

savefig(pl_ch, savedir * "lifetime_distributions_r_charting_time_shapes_recolor_lowr.png")

#----------------------------rp---------------------------------
data_dir_traj_root = "../../data/trajs_precomp_test6/"
data_dir_lt_root = "../../data/lifetimes_precomp_test6/"
data_dir_lt = data_dir_lt_root * "rp/"
data_dir_traj = data_dir_traj_root * "rp/"

files_trajs = readdir(data_dir_traj)
rs = vec(readdlm(data_dir_traj * files_trajs[1]))
rc = parse(Float64, files_trajs[1][end-5:end-4])

plots_topk = []
for (k_idx, k) in enumerate(ks)

    data_files = filter(f -> occursin("Δt1.", f) && occursin("k_$(k)_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc) r_p = ")

    push!(plots_topk, pl_distrib_x)
end

pl_topk = plot(plots_topk..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)
savefig(pl_topk, savedir * "lifetime_distributions_rp_rc_mismatch_change_rp_topk_lowr.png")

plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    data_files = filter(f -> occursin("Δt$(Δt_chart).txt", f) && occursin("k_10_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc) r_p = ")

    push!(plots_charting, pl_distrib_x)
end

pl_ch = plot(plots_charting..., layout=(3, 1), size=figsize, legend=:bottomleft, legendfontsize=legendfontsize, dpi=300)
savefig(pl_ch, savedir * "lifetime_distributions_rp_rc_mismatch_change_rp_chart_time_lowr.png")

#------------------------rc----------------------------

data_dir_traj_root = "../../data/trajs_precomp_test6/"
data_dir_lt_root = "../../data/lifetimes_precomp_test6/"
data_dir_lt = data_dir_lt_root * "rc/"
data_dir_traj = data_dir_traj_root * "rc/"

files_trajs = readdir(data_dir_traj)
rs = vec(readdlm(data_dir_traj * files_trajs[1]))
rp = parse(Float64, files_trajs[1][end-4:end-3])

plots_topk = []
for (k_idx, k) in enumerate(ks)

    data_files = filter(f -> occursin("Δt1.", f) && occursin("k_$(k)_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_p = %$(rp) r_c = ")

    push!(plots_topk, pl_distrib_x)
end

pl_topk = plot(plots_topk..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)
savefig(pl_topk, savedir * "lifetime_distributions_rp_rc_mismatch_change_rc_topk_lowr.png")


#keep k the same, change charting time
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    data_files = filter(f -> occursin("Δt$(Δt_chart).txt", f) && occursin("k_10_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_p = %$(rp) r_c = ")

    push!(plots_charting, pl_distrib_x)
end

pl_ch = plot(plots_charting..., layout=(3, 1), size=figsize, legend=:bottomleft, legendfontsize=legendfontsize, dpi=300)
savefig(pl_ch, savedir * "lifetime_distributions_rp_rc_mismatch_change_rc_charting time.png")

#=

#---------------------------alpha------------------------------

plots_topk = []
for (k_idx, k) in enumerate(ks)

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:none,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"\alpha = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#
    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end


pl_topk = plot(plots_topk..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)
savefig(pl_topk, savedir * "lifetime_distributions_alpha_memory_topk.png")


#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#READ trajectories
data_dir = data_dir_root * "alpha/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)[2:end]
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end


#------------params---------------

#gs4 = rs4 ./ α
#αs4 = rp ./ gs4

αs4 = [0.01, 0.0085, 0.0025, 0.00175]

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]



lognormal_intervals = [100:1:2000,100:1:2000,100:1:2000]
gaussian_intervals = [1:1:1000,100:1:1200,100:1:1200]

#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)
    push!(lifetimes_x_all_cht, lifetimes_x)
    push!(lifetimes_y_all_cht, lifetimes_y)
end

#keep k the same, change charting time
plots_charting = []
for (Δt_idx,Δt_chart) in enumerate(Δts_chart) 

    lifetimes_x = lifetimes_x_allk[Δt_idx]
    lifetimes_y = lifetimes_y_allk[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:none,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"\alpha = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts)*Δt_chart,std(lts.*Δt_chart)) 
        plot!(pl_distrib_y,gaussian_intervals[Δt_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end


pl_ch = plot(plots_charting..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)
savefig(pl_ch, savedir*"lifetime_distributions_alpha_memory_charting time.png")


#----------------------------change c parameter-------------------------
#increase c -> increase competition 

#READ trajectories
data_dir = data_dir_root * "/c/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end

cs = [1.5:0.3:2.5;]
lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx, k) in enumerate(ks)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, k, Δts_chart[1])
    push!(lifetimes_x_allk, lifetimes_x)
    push!(lifetimes_y_allk, lifetimes_y)
end

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")


    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=cs,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:none,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"c = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end


pl_topk = plot(plots_topk...,layout=(3,1),size=figsize,dpi=300)
savefig(pl_topk, savedir*"lifetime_distributions_competition_topk.png")

#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)
    push!(lifetimes_x_all_cht, lifetimes_x)
    push!(lifetimes_y_all_cht, lifetimes_y)
end

#keep k the same, change charting time
nbins=40
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_allk[Δt_idx]
    lifetimes_y = lifetimes_y_allk[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=cs, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:none,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"c = ")

    #plot normal distributions 
    for (i, lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts) * Δt_chart, std(lts .* Δt_chart))
        plot!(pl_distrib_y, gaussian_intervals[Δt_idx], gaussian, lc=colors_y[i], label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x, pl_distrib_y, layout=(1, 2), size=(800, 400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end


pl_ch = plot(plots_charting..., layout=(3, 1), size=figsize, dpi=300)
savefig(pl_ch, savedir*"lifetime_distributions_competition_charting time.png")
=#