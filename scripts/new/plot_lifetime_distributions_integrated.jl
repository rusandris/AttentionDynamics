using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")


#---------------------------lifetimes from integrated values------------------------------------ 

#read file names
#data_dir_traj_root = "../../data/trajs_precomp_test_c/"
#data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_traj_root = "../../data/trajs_precomp_alpha/"
data_dir_traj = data_dir_traj_root 
files = readdir(data_dir_traj)
files_trajs = files[findall(f -> occursin(".dat",f),files)] #select binary only

#TODO: save only x variables into  binary


#read params
#cs = vec(readdlm(data_dir_traj * files[1]))
#rs = vec(readdlm(data_dir_traj * files[1]))
αs = vec(readdlm(data_dir_traj * files[1]))
r = parse(Float64, files_trajs[1][10:11])
Δt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
T = parse(Float64, files_trajs[1][end-25:end-21])
#Δt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
#T = parse(Float64, files_trajs[1][end-25:end-21])
t = [0.0:Δt:T;]

plots_expbin = [plot(),plot()]
plots_linbin = [plot(), plot()]

colors = [:gray10, :orange,:darkgreen, :purple4,:lightblue,:darkblue]
markershapes = [:circle,:cross]

nbins = 30 #[25,10]
chart_times = [0.2,2.0] #[1.0,2.0]
data_shape = (Int(T/Δt)+1,300)
k = 10

for (i,f) in enumerate(files_trajs)
    @show f
    @show αs[i]

    #read trajectory
    traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))
    color = colors[i]

    for (j,chart_time) in enumerate(chart_times)
        @show chart_time

        #integrate activity segment to segment 
        int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt))
        lts = get_lifetimes_integ(int_act, k)

        #get exponential binning automatically adjusted for each
        bins_exp, counts_exp = autohist(lts; m=1, nbins=nbins)

        #get linear binning automatically adjusted for each
        h = fit(Histogram, lts; nbins=nbins)

        bins_lin = h.edges[1]
        counts_lin = h.weights

        #make plots 
        plot!(plots_expbin[j], log10.(bins_exp[1:end-1]), log10.(counts_exp), markershape=markershapes[j],
            mc=color, markerstrokecolor=color, 
            lc=color, 
            dpi=300,
            title=L"k=%$k", ms=4,
            label=L"\alpha = %$(αs[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")

        plot!(plots_linbin[j], log10.(bins_lin[1:end-1]), log10.(counts_lin), markershape=markershapes[j],
            mc=color, markerstrokecolor=color,
            lc=color,
            dpi=300,
            title=L"k=%$k", ms=4,
            label=L"\alpha = %$(αs[i]), Δt=%$(chart_time)",
            xlabel=L"\log(\tau)",
            ylabel=L"\log(count)")
    end
end
pl_expbin = plot(plots_expbin..., layout=(1, 2),size=(800,300), plot_title="expbin",
    legend=:bottomleft, left_margin=6Plots.mm)

pl_linbin = plot(plots_linbin..., layout=(1, 2), size=(800, 300), plot_title="linbin",
    legend=:bottomleft, left_margin=6Plots.mm)

pl_all = plot(pl_expbin, pl_linbin,layout=(2,1),size=(800,800))

plot!(pl_loglog[2],legend=:topright)

savedir = "../../paper_figs/lifetime_distributions/test_alpha/"
mkpath(savedir)
savefig(pl_all, savedir * "lifetime_distributions_k10_change_alpha_charting time_0.2_2_integrated_log_log_explin_bin_separated_properbins.png")

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

#=

ex,m,x_topk = get_inequality(traj,10,5)

pl_topk = plot(legend=false,size=(1000,300),dpi=300,ylims=(0.0,0.07),
    ylabel=L"x_i",xlabel=L"t",guidefontsize=14,left_margin=6Plots.mm,bottom_margin=6Plots.mm)

for i in topk_indxs
    plot!([ts[i]],x_topk[i]',st=:scatter,mc=colors')
    vline!([ts[i]],lc=:gray50,lw=1,la=0.3)
end

=#