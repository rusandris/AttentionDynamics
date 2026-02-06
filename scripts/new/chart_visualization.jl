using DelimitedFiles
using Plots, LaTeXStrings
cd(@__DIR__)
include("chart_functions.jl")


#---------------------------read trajectory files------------------------------------ 

#read file names
#data_dir_traj_root = "../../data/trajs_precomp_test_c/"
#data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_traj_root = "../../data/trajs_precomp_alpha/"
data_dir_traj = data_dir_traj_root
files = readdir(data_dir_traj)
files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only

#read an example traj

#read params
αs = vec(readdlm(data_dir_traj * files[1]))
r = parse(Float64, files_trajs[1][10:11])
dt = parse(Float64, files_trajs[1][end-6:end-4]) #sampling of trajectories
T = parse(Float64, files_trajs[1][end-25:end-21])
t = [0.0:dt:T;]
N = 300
data_shape = (Int(T / dt) + 1, N)

#----------------------------read file---------------------
file_example = files_trajs[1]
traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * file_example)), data_shape))

#decide timespan 
tspan = (1,20)
tspan_end_idx = Int(tspan[end]/dt)
tspan_start_idx = Int(tspan[1]/dt) + 1 
tidxs = tspan_start_idx:tspan_end_idx

#compute integrated variables
int_act, time_segments = integrate_activity(traj, t, Int(Δt/dt));
int_act_matrix = stack(int_act)'
t_int = [time_segments[i][end] for i in 1:length(time_segments)]

tspan_end_intidx = Int(tspan[end] / Δt)
tspan_start_intidx = Int(tspan[1] / Δt) 
tint_idxs = tspan_start_intidx:tspan_end_intidx



#------------------------------visualize charting process-----------------------------------

#find out who is active in that timespan
k = 5
topk_idxs, topk_vals, L_minmax, L_mean = chart_stats(int_act_matrix[tint_idxs, :], k, 1)

#count number of colors needed
idxs_active = unique(stack(topk_idxs))
nr_colors = length(idxs_active)
#colors = range(HSV(0, 1, 1), stop=HSV(-300, 1, 1), length=26)
#range(HSV(colorant"red"), stop=HSV(colorant"yellow"), length=nr_colors)
colors = distinguishable_colors(nr_colors, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
idx_to_color = Dict(zip(idxs_active,colors))

#plot trajectory itself
pl_ts = plot(t[tidxs], traj[tidxs, idxs_active], lc=colors', lw=3,
    ylabel=L"L_i(t)", xlabel="", legend=false, xformatter=:none,
    tickfontsize=6, guidefontsize=18, left_margin=5Plots.mm)

#add vertical lines 
Δt = 1.0
n_Δt = Int(Δt / dt)
vline!(pl_ts, t[tidxs[1]:n_Δt:tidxs[end]], ls=:dash, lc=:gray50, la=0.5)

#plot inactive with gray separately
idxs_inactive = []
for i in 1:N
    i ∈ idxs_active && continue
    push!(idxs_inactive,i)
end

plot!(pl_ts,t[tidxs], traj[tidxs, idxs_inactive], lc=:gray50, lw=1,
    ylabel=L"L_i(t)", xlabel="", legend=false, xformatter=:none,
    tickfontsize=10, guidefontsize=18, left_margin=5Plots.mm)

#plot integrated variables
pl_int_act = plot(t_int[tint_idxs], int_act_matrix[tint_idxs, idxs_active],
    lc=colors', mc=colors', lw=0.5,markershape=:circle,ms=6,
    markerstrokewidth=0.0,
    ylabel = L"\tilde{L}_i(t,\Delta t)",xlabel="",
    legend=false, tickfontsize=6, guidefontsize=15,
    left_margin=5Plots.mm, xformatter=:none)

#ylabel = L"\int_{Δt_i}^{Δt_{i+1}} L_i(t)dt"
#add vertical lines 
Δt = 1.0
n_Δt = Int(Δt / dt)
vline!(pl_int_act, tint_idxs, ls=:dash, lc=:gray50, la=0.5)

#plot inactive with gray separately
plot!(t_int[tint_idxs], int_act_matrix[tint_idxs, idxs_inactive],
    lc=:gray50, mc=:gray50, lw=1,
    legend=false, tickfontsize=10, guidefontsize=15,
    left_margin=5Plots.mm, xformatter=:none)

#plot charts and ranks
pl_topk = plot(legend=false,dpi=300,
    ylabel="ranks",xlabel=L"t",guidefontsize=18,
    left_margin=5Plots.mm,bottom_margin=6Plots.mm,
    yflip=true, yticks=[1:k;], ylims=(0, k + 1),
    tickfontsize=15)

cons_occ = consecutive_occurrences(topk_idxs)

for i in 1:length(cons_occ)
    color = idx_to_color[idxs_active[i]]
    plot!(pl_topk, t_int[tint_idxs], cons_occ[i],
        markerstrokecolor=color,
        mc=color, ms=6,lc=color,
        markershape=:circle,lw=0.5)
end

#=
for (i,t) in enumerate(t_int[tint_idxs]) 
    chart = topk_idxs[i]
    
    colors_chart = [idx_to_color[i] for i in chart]
    plot!(pl_topk, [t], [1:k;]', st=:scatter, 
        markerstrokecolor=colors_chart', 
        mc=colors_chart', ms=6)
    
end
=#
plot!(pl_ts,ylims=(0.0,0.03),yticks=[0.0,0.015,0.03])
plot!(pl_int_act, ylims=(0.0, 0.03), yticks=[0.0, 0.015, 0.03])

pl = plot(pl_ts, pl_int_act, pl_topk, layout=(3, 1), size=(800, 600),dpi=300)
savefig(pl,"../../paper_figs/chart_explanation_example.png")