using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")


#---------------------------lifetimes from integrated values------------------------------------ 

#read file names
#data_dir_traj_root = "../../data/trajs_precomp_test_c/"
data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_traj = data_dir_traj_root * "rp=rc/"
files = readdir(data_dir_traj)
files_trajs = files[findall(f -> occursin(".dat",f),files)][[4,1,2,3]] #select binary only

#TODO: save only x variables into  binary


#read params
#cs = vec(readdlm(data_dir_traj * files[1]))
rs = vec(readdlm(data_dir_traj * files[1]))
#r = parse(Float64, files_trajs[1][10:12])
#Δt = parse(Float64, files_trajs[1][end-7:end-4]) #sampling of trajectories
#T = parse(Float64, files_trajs[1][end-26:end-22])
Δt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
T = parse(Float64, files_trajs[1][end-25:end-21])
t = [0.0:Δt:T;]

pl_loglog = plot();
pl_linlin = plot();
colors = [[:gray10, :gray40],[:red, :orange],[:lightgreen,:darkgreen],
    [:skyblue,:darkblue]]

nbins = 20
chart_times = [0.5,3.0] #[1.0,2.0]
data_shape = (Int(T/Δt)+1,300)
k = 20

for (i,f) in enumerate(files_trajs)
    @show f
    @show rs[i]

    #read trajectory
    traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))
    colorset = colors[i]

    for (j,chart_time) in enumerate(chart_times)
        @show chart_time

        #integrate activity segment to segment 
        int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt))
        lts = get_lifetimes_integ(int_act, k)
        @show length(lts)
        @show extrema(lts)

        plot_distrib!(pl_loglog, lts, 1.0;
            color=colorset[j], nbins=nbins, mode=:none,
            scale=:log,binning=:exp,st=:path,
            title=L"k=%$k", ms=4,
            label_par=L"r = %$(rs[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot_distrib!(pl_linlin, lts, 1.0;
            color=colorset[j], nbins=nbins, mode=:none,
            scale=:lin, binning=:lin, st=:path,
            title=L"k=%$k", ms=4,
            label_par=L"r = %$(rs[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel="count")

    end
end

savedir = "../../paper_figs/lifetime_distributions/test8_origpar/"
mkpath(savedir)
savefig(pl_loglog, savedir * "lifetime_distributions_k10_change_rp=rc_charting time_0.5_3_integrated_log_log_exp_bin.png")

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
