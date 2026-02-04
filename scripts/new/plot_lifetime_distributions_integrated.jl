using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")


#---------------------------lifetimes from integrated values------------------------------------ 

#read file names
#data_dir_traj_root = "../../data/trajs_precomp_test_c/"
#data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_traj_root = "../../data/trajs_precomp_test_c/"
data_dir_traj = data_dir_traj_root * "c_low/"
files = readdir(data_dir_traj)
files_trajs = files[findall(f -> occursin(".dat",f),files)] #select binary only

#TODO: save only x variables into  binary


#read params
cs = vec(readdlm(data_dir_traj * files[1]))
#rs = vec(readdlm(data_dir_traj * files[1]))
r = parse(Float64, files_trajs[1][10:11])
Δt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
T = parse(Float64, files_trajs[1][end-25:end-21])
#Δt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
#T = parse(Float64, files_trajs[1][end-25:end-21])
t = [0.0:Δt:T;]

pl_loglog = plot();
pl_linlin = plot();
colors = [[:gray10, :gray40],[:red, :orange],[:lightgreen,:darkgreen],
    [:skyblue, :darkblue], [:maroon3, :maroon], [:plum4, :purple4]]

nbins = [25,10]
chart_times = [0.5,3.0] #[1.0,2.0]
data_shape = (Int(T/Δt)+1,300)
k = 10

#predefine bin edges for all datasets
lt_min = 1
lt_max = 200
b = lt_max^(1 / nbins[1]) #base 
bin_edges_exp =  b .^ [0:1:nbins[1];]
bin_edges_lin = range(lt_min, lt_max; length=nbins[1])

for (i,f) in enumerate(files_trajs)
    @show f
    @show cs[i]

    #read trajectory
    traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))
    colorset = colors[i]

    for (j,chart_time) in enumerate(chart_times)
        @show chart_time
        nbin = nbins[j]

        #integrate activity segment to segment 
        int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt))
        lts = get_lifetimes_integ(int_act, k)
        @show length(lts)
        @show extrema(lts)

        #get counts and bins
        counts, bins = get_distrib(τs, bin_edges_exp; scale=:log)

        #make plots 
        plot!(pl_loglog, counts, bins, markershape=:circle,
            mc=color, markerstrokecolor=color, 
            lc=color, 
            label=label_par, dpi=300,
            title=L"k=%$k", ms=4,
            label_par=L"c = %$(cs[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        counts, bins = get_distrib(τs, bin_edges_lin; scale=:lin)

        plot!(pl_linlin, counts, bins, markershape=:circle,
            mc=color, markerstrokecolor=color, 
            lc=color, 
            label=label_par, dpi=300,
            title=L"k=%$k", ms=4,
            label_par=L"c = %$(cs[i]), Δt=%$(chart_time)",
            xlabel=L"\tau",
            ylabel=L"count")

    end
end

savedir = "../../paper_figs/lifetime_distributions/test_c_small/"
mkpath(savedir)
savefig(pl_loglog, savedir * "lifetime_distributions_k10_change_c_charting time_0.5_3_integrated_log_log_exp_bin.png")

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
