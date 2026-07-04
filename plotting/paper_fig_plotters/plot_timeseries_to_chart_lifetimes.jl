using Plots, LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

#-------------------------------N=10(LY)------------------------------------

datadir_traj_root = "../../data/trajs_precomp_orig_N10/"
data_dir_traj = datadir_traj_root * "traj_N10_B1B2B3B4C_Δt0.01_T5000_Ttr1000_examples/"

N, Δt_traj, T = readdlm(data_dir_traj * "traj_spec_params.txt")
Ttr = 1000
Δt_traj = 1e-2

files_trajdata = readdir(data_dir_traj)[2:end]
nr_ics = 1
files_trajdata_examples = files_trajdata[1:nr_ics:end]
data_shape = (Int(T * 1/Δt_traj + 1), Int(N))
t_plot = 20
ts = [Ttr:Δt_traj:Ttr+T;]

titles = [L"B_1",L"B_2",L"B_3",L"B_4",L"C"]
title_colors = [:cyan :lime :orange :purple]
guidefontsize = 30
legendfontsize = 12
tickfontsize = 14
titlefontsize = 24
annot_fontsize = 30

colors = distinguishable_colors(300,
    [RGB(1, 1, 1), RGB(0, 0, 0), RGB(0.8, 1.0, 0.8), RGB(0.6274, 1.0, 1.0), RGB(1.0, 1.0, 0.0)],
    dropseed=true, lchoices=range(0, stop=80, length=15),)

plots_ts_N10 = []
yticks = [[0.0, 0.05, 0.1], [0.0, 0.1, 0.2], [0.0, 0.1], [0.0, 0.1, 0.2, 0.3]]


for (i, f) in enumerate(files_trajdata_examples[1:end-1])
    @show f
    traj = reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape)

    colors_timeseries = colors[1:10] 
    interval = 1:Int(t_plot * 1/Δt_traj) 
    pl = plot(ts[interval], traj[interval,:], legend=false, ylabel=L"L_i",
        guidefontsize=guidefontsize-4, legendfontsize=legendfontsize,
        tickfontsize=tickfontsize, titlefontsize=titlefontsize,lw=2,
        lc=colors_timeseries', title=titles[i], titlefontcolor=title_colors[i],
        yticks=yticks[i])

    #annotate!(pl,Ttr+1, 0.3, text(titles[i], title_colors[i], :center, annot_fontsize))
    
    push!(plots_ts_N10,pl)
end

#adjust plots 
plot!(plots_ts_N10[end], xlabel=L"t_L", guidefontsize=guidefontsize)
#plot!(plots_ts_N10[1],title = L"N = 10")

pl_N10 = plot(plots_ts_N10..., layout=(4, 1), size=(1000, 1200), left_margin=5Plots.mm,ylims=(0, 0.32))


#-------------------------------N=300------------------------------------

datadir_traj_root = "../../data/trajs_precomp_orig/"
#data_dir_traj = datadir_traj_root * "horizontal_and_curve_article_test_Δt0.1_T5000_Ttr500/"
data_dir_traj = datadir_traj_root * "trajs_examples_B1B2B3B4C_Δt0.02_T1000_Ttr2000/"

N, Δt_traj, T = readdlm(data_dir_traj * "traj_spec_params.txt")
Ttr = 2000
N = Int(N)
#Δt_traj = 1e-1

files_trajdata = readdir(data_dir_traj)[2:end]
nr_ics = 1
files_trajdata_examples = files_trajdata[1:nr_ics:end]
data_shape = (Int(T * 1 / Δt_traj + 1), Int(N))
t_plot = 20
ts = [Ttr:Δt_traj:T+Ttr;]

plots_ts_N300 = []
yticks = [[0.0, 0.015, 0.03], [0.0, 0.1, 0.2], [0.0, 0.1], [0.0, 0.05, 0.1]]


for (i, f) in enumerate(files_trajdata_examples[1:end-1])
    @show f
    traj = reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape)

    interval = 1:Int(t_plot * 1 / Δt_traj)

    #select actives 
    traj_segm = traj[interval, :]
    idxs_active = [1:N;][[maximum(traj_segm[:,i]) > 1e-2 for i in 1:N]]
    colors_active = colors[1:length(idxs_active)]
    colors_timeseries = fill(RGB(0.9),N)
    colors_timeseries[idxs_active] .= colors_active
    las = fill(0.1,N)
    las[idxs_active] .= 1.0
    lws = fill(0.5, N)
    lws[idxs_active] .= 2

    pl = plot(ts[interval], traj[interval, :], legend=false, ylabel="",
        guidefontsize=guidefontsize-4, legendfontsize=legendfontsize,
        tickfontsize=tickfontsize, titlefontsize=titlefontsize, lw=lws',
        lc=colors_timeseries', title=titles[i], titlefontcolor=title_colors[i],
        yticks=yticks[i])

    #annot_x = i > 1 ? 0.2 : 0.04
    #annotate!(Ttr + 1, annot_x, text(titles[i], title_colors[i], :center, annot_fontsize))

    push!(plots_ts_N300, pl)
end

#adjust plots 
plot!(plots_ts_N300[end], xlabel=L"t_L", guidefontsize=guidefontsize)
pl_N300 = plot(plots_ts_N300..., layout=(4, 1), size=(1000, 1200), ylims=(0, 0.22))

#add title (annotation)
annotate!(0.5, 0.0, text(L"N = 300", :gray10, :center, 25))

plot!(pl_N300[1], ylims=(0.0, 0.05),title = L"N=300")
plot!(pl_N300,yticks=[0,0.1,0.2])
plot!(pl_N300[1], yticks=[0, 0.02, 0.04])

pl_log300 = plot(pl_N300, yscale=:log10, ylabel="", ylims=(1e-10, 1e2),
    yticks=10.0 .^ [2:-4:-10;])

pl_linlog300 = plot(pl_N300, pl_log300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm, dpi=300)

savefig(pl_linlog300, "../../paper_figs/timeseries_examplesB1B2B3B4_N300_T1000_Ttr2000_Δt0.02_side-by-side_linlog.png")

#N10 plots

pl = plot(pl_N10, pl_N300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm, dpi=300)

pl_logN300 = plot(pl_N300, yscale=:log10, ylabel="", ylims=(1e-10, 1e2),
    yticks=10.0 .^ [2:-4:-10;])

pl_linlog300 = plot(pl_N300, pl_logN300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm, dpi=300)

pl_logN10 = plot(pl_N10, yscale=:log10, ylabel="", ylims=(1e-10, 1e2),
    yticks=10.0 .^ [2:-4:-10;])

pl_linlog10 = plot(pl_N10, pl_logN10, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm, dpi=300)

pl_linN10_N300 = plot(pl_N10, pl_N300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm, dpi=300)

pl_logN10_N300 = plot(pl_logN10, pl_logN300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm,
    left_margin=10Plots.mm,  dpi=300)

savefig(pl_linN10_N300, "../../paper_figs/timeseries_examplesB1B2B3B4_N10_N300_T1000_Ttr2000_Δt0.02_side-by-side_lin.png")
savefig(pl_linN10_N300, "../../paper_figs/timeseries_examplesB1B2B3B4_N10_T1000_Ttr2000_Δt0.02_side-by-side_lin.png")
savefig(pl_logN10_N300, "../../paper_figs/timeseries_examplesB1B2B3B4_N10_T1000_Ttr2000_Δt0.02_side-by-side_log.png")


#savefig(pl,"../../paper_figs/timeseries_timescales_N10_T_Ttr_N300_N300_T1000_Ttr2000_Δt0.02_side-by-side.png")
#savefig(pl_log, "../../paper_figs/timeseries_timescales_N10_T_Ttr_N300_N300_T1000_Ttr2000_Δt0.02_side-by-side_loglog.png")