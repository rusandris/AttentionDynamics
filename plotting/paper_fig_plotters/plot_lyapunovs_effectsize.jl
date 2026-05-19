using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#path to lifetimes
data_dir_lt_root = "../../data/lyapunovs_effect_size/"
#data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws461_T_WLC1000_T_fp10_Ttr_WLC10000_Ttr_fp10_dtmax_0.001/"
data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws_largerN_T_WLC1000_T_fp10_Ttr_WLC1000_Ttr_fp10_dtmax_0.001/"
data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws_largerN_gnormalized_phasediag_sec1_T_WLC1000_T_fp10_Ttr_WLC1000_Ttr_fp10_dtmax_0.001/"

files = readdir(data_dir_lt)
files_lyaps = readdir(data_dir_lt)[1:end-1][[2, 3, 4, 1]] #[[4,5,6,1,2,3]]
Ns = Int.(readdlm(data_dir_lt * "system_sizes.txt"))

#plot all on one fig 

pl = plot(ylabel=L"\lambda", xlabel=L"w",guidefontsize=20,tickfontsize=15,
    legendfontsize = 12,legend=:topright,legendcolumns=3,ylims=(-0.5,4.5),
    size=(1200, 400), bottom_margin=8Plots.mm, right_margin=8Plots.mm, 
    left_margin=8Plots.mm,dpi=300,xlims = (0.0,1.0),xticks=[0.2:0.2:1.0;])
colors = distinguishable_colors(length(Ns)+3)[3:end]

wc = 1/(2000/299) #1/4
vline!([wc], ls=:dash, lc=:gray10,lw=3, label=false)

for (i,file) in enumerate(files_lyaps)
    N = Ns[i]
    #read lyap exps 
    lyap_data = readdlm(data_dir_lt * file)
    ws = lyap_data[:,1]
    λs =  lyap_data[:,2]
    plot!(pl, ws, λs, lc = colors[i], lw=2, label=L"N = %$(N)")
    vline!([(N - 1) / N], ls=:dash, lc=colors[i], label="")

end
savefig(pl,"../../paper_figs/lyapunovs_system_size_effect_chaotic_width_Ttr1e3_tscale6.688_crop.png")

#example timeseries 
using OrdinaryDiffEq
using DynamicalSystems
N_test = 300
u0_test = rand(2N_test)
p_test = [2000 / 299 * (N_test-1), 0.25, N_test]
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)

ds_test = create_ds(u0_test, p_test; diffeq=diffeq)
traj1, t = trajectory(ds_test, 100, u0_test; Ttr=1000, Δt=1e-3)

pl1 = plot(t[1:1000], Matrix(traj1)[1:1000, 1:2:end], legend=false,
    xlabel=L"t", ylabel=L"x_i",title=L"N=%$(N_test),w=0.25",guidefontsize=25)

diffeq2 = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
p_test2 = [2000 / 299 * (N_test - 1), 0.4, N_test]
ds_test2 = create_ds(u0_test, p_test2; diffeq=diffeq2)
traj2, t = trajectory(ds_test2, 100, u0_test; Ttr=2000, Δt=1e-3)
pl2 = plot(t[1:1000], Matrix(traj2)[1:1000, 1:2:end], legend=false,
    xlabel=L"t", ylabel=L"x_i", title=L"N=%$(N_test),w=0.4", guidefontsize=25)

pl = plot(pl1, pl2, layout=(1, 2), size=(800, 400), left_margin=6Plots.mm, bottom_margin=6Plots.mm,
    right_margin=6Plots.mm)
#savefig(pl,"../../paper_figs/timeseries_examples_N100_lyapunovs_system_size_effect_chaotic_width_Ttr1e3.png")

#L system
u0_test = rand(2N_test) .* 1e-2
p_test = [10, 10, 1.0, 2.4, 0.005, 300]
ds_L = create_original_ds(u0_test, p_test; diffeq=diffeq2)
traj_L, t = trajectory(ds_L, 100, u0_test; Ttr=100, Δt=1e-3)
pl_L = plot(t, Matrix(traj_L)[:, 1:2:end], legend=false,
    xlabel=L"t", ylabel=L"L_i", title=L"N=%$(N_test),L", guidefontsize=25)

#x system
g, w = rescaled_params(p_test...)
ds_x = create_ds(u0_test, [g,w,300]; diffeq=diffeq2)

#init cond
η_L,η_Y = rescaling_factors(p_test...)
u0_test[1:2:end] .*= 1 / η_L
u0_test[2:2:end] .*= 1 / η_Y
traj_x, t = trajectory(ds_x, 10, rand(2N_test); Ttr=0, Δt=1e-3)
