using AttentionDynamics
using Revise
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

data_dir_traj = "../../data/trajs_precomp_orig/"

#---------------------------read trajectories filenames------------------------------------ 
dir_content = readdir(data_dir_traj)
files = dir_content[findall(f -> occursin(".txt", f) || occursin(".dat", f), dir_content)]

files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only

#files corresponding to a certain param set
par_names = ["B1", "B2", "C", "B3", "B4"]
file_sets = [files_trajs[findall(f -> occursin(par_names[i], f), files_trajs)] for i in 1:5]

#lifetime settings 
k = 10
chart_times = [0.1, 0.5, 1.0] #2.0]
Δt_traj = 0.1
data_shape = (50010,300) 
for fileset in file_sets
    f = fileset[1] #test first of every fileset
    traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))
    n_charting = Int(chart_times[2] / Δt_traj)
    charts, _ = chart_stats(traj, k, n_charting; threshold=1e-3)
    p_one = probability_first_entry(charts)
    @show p_one
end