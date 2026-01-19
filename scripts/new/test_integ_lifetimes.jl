include("chart_functions.jl")
cd(@__DIR__)


nbins = 50
figsize = (600, 800)
legendfontsize = 10
data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_lt_root = "../../data/lifetimes_precomp_test_integ_origpar/"


#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt
N=300
Δt = 1e-1
K = 1.0;
c = 2.4;
α = 0.005;

#READ trajectories
data_dir = data_dir_traj_root * "rp/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)[2:end][[4,1,2,3]]
trajs = Matrix{Float64}[]

for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end
#read params 
#cs = vec(readdlm(data_dir * files_all[1]))
rs = vec(readdlm(data_dir * files_all[1]))
#r = parse(Float64,traj_files[1][16:18])
#r = parse(Float64, traj_files[1][12:14])

#get lifetimes with different k, charting times 
ks = [10, 50, 100]
ns_chart = [10, 20, 50]
ts = collect(500:0.1:10500)

#precompute lifetimes 
savedir = data_dir_lt_root * "rp=rc/"
mkpath(savedir)

for (n_idx, n) in enumerate(ns_chart)

    for i in 1:length(trajs)
        traj = trajs[i][:,1:2:2N]
        L_ints,time_segments = integrate_activity(traj, ts, n)
        lts = get_lifetimes_integ(L_ints, ks[1])
        @show length(lts)
        writedlm(savedir * "lifetimes_integ_rp=rc=$(rs[i])_c=$(c)_k_$(ks[1])_Δt_int$(n).txt", lts)

    end

end

#-------------------test if it can be read and plotted-------------------------

using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")
#-------------------------plot settings------------------------------

#get lifetimes with different k, charting times
savedir = "../../paper_figs/lifetime_distributions/test_integ_origpar/"
mkpath(savedir)

nbins = 50
figsize = (400, 800)
legendfontsize = 8
Δt = 1e-1 #sampling rate of trajs 

#------------------------------plot test------------------------------
data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_lt_root = "../../data/lifetimes_precomp_test_integ_origpar/"


data_dir_lt = data_dir_lt_root * "rp=rc/"
data_dir_traj = data_dir_traj_root * "rp=rc/"
files_trajs = readdir(data_dir_traj)

rs = vec(readdlm(data_dir_traj * files_trajs[1]))

plots_int = []
for (n_idx, n) in enumerate(ns_chart)

    data_files = filter(f -> occursin("int$n", f) && occursin("k_10_", f), readdir(data_dir_lt))
    @show length(data_files)
    lifetimes_x = []
    for f in data_files
        lts = readdlm(data_dir_lt * f)
        push!(lifetimes_x, lts)
    end

    @show length(lifetimes_x)
    colors_x = distinguishable_colors(length(lifetimes_x))

    pl_distrib_x = plot_distribs(lifetimes_x, n;
        ps=rs, colors=colors_x, nbins=nbins, mode=:none,
        scale=:log,
        title=L"k=10" * L",Δt=%$(n*Δt)",
        label_par=L"rp = = rc =")

    push!(plots_int, pl_distrib_x)
end

pl_int = plot(plots_int..., layout=(3, 1), size=figsize,
    legend=:topright,legendfontsize=legendfontsize, dpi=300)

savefig(pl_int, savedir * "lifetime_distributions_integ_rp_rc_orig_params.png")


