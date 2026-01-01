using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics
using OrdinaryDiffEq
using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

#------------------IDEA---------------------
#test multidim OP after critical point 
#what happens to the symbolic dynamics of w is changed?
#PE,WPE
#distribution / average lifetime of in the top k list 
#changing a clearly interpretable parameter: α,c 
#corresponds to moving along a curve in g-w space
#which parameter change leads to driving the system more towards criticality?
#how does properties of the system change when it gets near criticality?


#TODO:
#try exponential/statistically significant binning 
#change charting time (3plot) + write time sampling on axis
#charting time 1/10th of average peak width
#add fits
#show inequality on fig

include("chart_functions.jl")

gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2))

nbins = 50
figsize = (600,800)
legendfontsize=12
savedir = "../../paper_figs/lifetime_distributions/"
mkpath(savedir)

#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

rs4 = [6.0:1.0:9.0;]
Δt = 1e-1
rc = 9.0;
rp = 9.0;
K = 1.0;
c = 2.4;
α = 0.005;

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rp=rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end


#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [1,10,30]

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

#plot lifetimes 

#colors_x = palette([:purple, :red], 4)
#colors_y = palette([:green, :blue], 4)
colors_x = distinguishable_colors(4)
colors_y = distinguishable_colors(4)

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

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

savefig(pl_topk,savedir*"lifetime_distributions_r_topk_shapes_recolor_lowr.png")


#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx,Δt) in enumerate(Δts_chart) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[1],Δt)
    push!(lifetimes_x_all_cht,lifetimes_x)
    push!(lifetimes_y_all_cht,lifetimes_y)
end


#keep k the same, change charting time
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_all_cht[Δt_idx]
    lifetimes_y = lifetimes_y_all_cht[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=rs4, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p=r_c = ")

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

pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), legendfontsize=legendfontsize, dpi=300)

savefig(pl_ch, savedir*"lifetime_distributions_r_charting_time_shapes_recolor_lowr.png")


# compare with rp \neq rc (mismatch between production and consumption rate)

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rp/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end

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

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=rs4,colors=colors_x,nbins=nbins,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins,mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",scale=:lin,
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

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
savefig(pl_topk, savedir*"lifetime_distributions_rp_rc_mismatch_change_rp_topk_lowr.png")

#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)
    push!(lifetimes_x_all_cht, lifetimes_x)
    push!(lifetimes_y_all_cht, lifetimes_y)
end

plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_all_cht[Δt_idx]
    lifetimes_y = lifetimes_y_all_cht[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=rs4, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p=r_c = ")

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


pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), legendfontsize=legendfontsize, dpi=300)
savefig(pl_ch, savedir * "lifetime_distributions_rp_rc_mismatch_change_rp_chart_time_lowr.png")


#WATCH OUT FOR NUMERICAL METHOD ERRORS
# compare with rp \neq rc (mismatch between production and consumption rate)
#ONLY CHANGE rc
#precompute trajectories 


#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end


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

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]


    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p = %$(rp), r_c = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=    
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x[1:2]) 
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
savefig(pl_topk, savedir*"lifetime_distributions_rp_rc_mismatch_change_rc_topk_lowr.png")

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
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")


    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p = %$(rp), r_c = ")

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
savefig(pl_ch, savedir*"lifetime_distributions_rp_rc_mismatch_change_rc_charting time.png")



#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/alpha/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
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

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
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
    push!(plots_topk,pl_distrib_x)
end


pl_topk = plot(plots_topk..., layout=(3, 1), size=figsize, legendfontsize=legendfontsize, dpi=300)
savefig(pl_topk, savedir*"lifetime_distributions_alpha_memory_topk.png")

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
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
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
data_dir = "../../data/trajs_precomp_lifetimes/c/"
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
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")


    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=cs,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
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
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=cs, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
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


pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), dpi=300)
savefig(pl_ch, savedir*"lifetime_distributions_competition_charting time.png")
