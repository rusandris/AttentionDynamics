using Plots, LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

#-------------------------------N=10(rescaled)------------------------------------

datadir_traj_root = "../../data/trajs_precomp_rescaled/"
data_dir_traj = datadir_traj_root * "wg_horizontal_and_curve_B1B2B3B4C_Δt0.01_T20000_Ttr500/"

N, Δt_traj, T = readdlm(data_dir_traj * "traj_spec_params.txt")
Ttr = 500
Δt_traj = 1e-2

files_trajdata = readdir(data_dir_traj)[3:end]
nr_ics = 10
files_trajdata_examples = files_trajdata[1:nr_ics:end]
data_shape = (Int(T * 1/Δt_traj + 1), Int(N))
t_plot = 10
ts = [Ttr:Δt_traj:T;]

titles = [L"B_1",L"B_2",L"B_3",L"B_4",L"C"]
guidefontsize = 30
legendfontsize = 12
tickfontsize = 14
titlefontsize = 24

plots_N10 = []


for (i, f) in enumerate(files_trajdata_examples)
    @show f
    traj = reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape)

    interval = 1:Int(t_plot * 1/Δt_traj) 
    pl = plot(ts[interval], traj[interval,:], legend=false, ylabel=L"x_i",
        guidefontsize=guidefontsize, legendfontsize=legendfontsize,
        tickfontsize=tickfontsize, titlefontsize=titlefontsize,lw=2,
        title = titles[i])
    
    push!(plots_N10,pl)
end

#adjust plots 
plot!(plots_N10[end], xlabel=L"t", guidefontsize=guidefontsize)

pl_N10 = plot(plots_N10..., layout=(5, 1), size=(1000, 1200), left_margin=5Plots.mm,
    plot_title=L"N=10 \textnormal{ (rescaled)}", plot_titlefontsize=25)


#-------------------------------N=300 (rescaled)------------------------------------

datadir_traj_root = "../../data/trajs_precomp_orig/"
#data_dir_traj = datadir_traj_root * "horizontal_and_curve_article_test_Δt0.1_T5000_Ttr500/"
data_dir_traj = datadir_traj_root * "trajs_for_chart_lifetime_examples_002/"

N, Δt_traj, T = readdlm(data_dir_traj * "traj_spec_params.txt")
Ttr = 500
#Δt_traj = 1e-1

files_trajdata = readdir(data_dir_traj)[2:end]
nr_ics = 1
files_trajdata_examples = files_trajdata[1:nr_ics:end]
data_shape = (Int(T * 1 / Δt_traj + 1), Int(N))
t_plot = 10
ts = [Ttr:Δt_traj:T+Ttr;]


plots_N300 = []

for (i, f) in enumerate(files_trajdata_examples)
    @show f
    traj = reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape)

    interval = 1:Int(t_plot * 1 / Δt_traj)
    pl = plot(ts[interval], traj[interval, :], legend=false, ylabel=L"x_i",
        guidefontsize=guidefontsize, legendfontsize=legendfontsize,
        tickfontsize=tickfontsize, titlefontsize=titlefontsize, lw=2,
        title=titles[i])

    push!(plots_N300, pl)
end

#adjust plots 
plot!(plots_N300[end], xlabel=L"t", guidefontsize=guidefontsize)

pl_N300 = plot(plots_N300..., layout=(5, 1), size=(1000, 1200),
    plot_title = L"N=300",plot_titlefontsize=25)

pl = plot(pl_N10, pl_N300, layout=(1, 2), size=(2000, 1200), bottom_margin=5Plots.mm, 
    left_margin=10Plots.mm,dpi=300)
savefig(pl,"../../paper_figs/timeseries_timescales_N10_N300_N300T500_Δt0.02_side-by-side.png")