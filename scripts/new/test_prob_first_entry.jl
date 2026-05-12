using AttentionDynamics
using Revise
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

data_dir_traj = "../../data/trajs_precomp_orig/horizontal_and_curve_article_test_Δt0.1_T5000_Ttr500/"

#---------------------------read trajectories filenames------------------------------------ 
dir_content = readdir(data_dir_traj)
files = dir_content[findall(f -> occursin(".txt", f) || occursin(".dat", f), dir_content)]

files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only

#files corresponding to a certain param set
par_names = ["B1", "B2", "C", "B3", "B4"]
file_sets = [files_trajs[findall(f -> occursin(par_names[i], f), files_trajs)] for i in 1:5]

#lifetime settings 
k = 5
chart_times = [0.1, 0.5, 1.0,2.0,10.0] #2.0]
Δt_traj = 0.1
data_shape = (50001,300)
t = [0.0:Δt_traj:5000;]
chart_time = chart_times[4]
n_charting = Int(chart_time / Δt_traj)
probs_one = []
for i in 1:length(par_names)
	@show par_names[i]
    f = file_sets[i][1] #test first of every fileset
    @show f
    
    #read traj
    traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))
    #integrated activity

    int_act, t_segments = integrate_activity(traj,t,n_charting)
    charts, _ = chart_stats(stack(int_act)', k, 1)
    p_one = probability_first_entry(charts)
    @show p_one
    push!(probs_one,p_one)
end

plot(bar([1:5],probs_one),xticks=(1:5,par_names),ylims=(0,1),
	ylabel=L"P_{one}",title=L"k=%$(k),\Delta t = %$(chart_time)")
