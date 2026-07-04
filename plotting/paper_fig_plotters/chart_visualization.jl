using DelimitedFiles
using Plots, LaTeXStrings
using AttentionDynamics
cd(@__DIR__)


#---------------------------read trajectory files------------------------------------ 

#read file names
#data_dir_traj_root = "../../data/trajs_precomp_test_c/"
#data_dir_traj_root = "../../data/trajs_precomp_test8_origpar/"
data_dir_traj_root = "../../data/trajs_precomp_orig_N300/"
data_dir_traj = data_dir_traj_root * "trajs_examples_B1B2B3B4C_Δt0.02_T1000_Ttr2000/"
files = readdir(data_dir_traj)
files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only

#read an example traj

#read params
N,Δt,T = vec(readdlm(data_dir_traj * files[1]))
t = [0.0:Δt:T;]
N = 300
data_shape = (Int(T / Δt) + 1, N)

#----------------------------read file---------------------
file_example = files_trajs[2]
traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * file_example)), data_shape))

#decide timespan 
tspan = (10,20)
tspan_end_idx = Int(tspan[end]/Δt)
tspan_start_idx = Int(tspan[1]/Δt) + 1 
tidxs = tspan_start_idx:tspan_end_idx

#compute integrated variables
Δt_chart = 1.0
int_act, time_segments = integrate_activity(traj, t, Int(Δt_chart/Δt));
int_act_matrix = stack(int_act)'
t_int = [time_segments[i][end] for i in 1:length(time_segments)]

tspan_end_intidx = Int(tspan[end] / Δt_chart)
tspan_start_intidx = Int(tspan[1] / Δt_chart) 
tint_idxs = tspan_start_intidx:tspan_end_intidx

#------------------------------visualize charting process-----------------------------------

#find out who is active in that timespan
k = 50
topk_idxs, _ = chart_stats(int_act_matrix[tint_idxs, :], k, 1)

#count number of colors needed
idxs_active = unique(stack(topk_idxs))
nr_colors = length(idxs_active)
#colors = range(HSV(0, 1, 1), stop=HSV(-300, 1, 1), length=26)
#range(HSV(colorant"red"), stop=HSV(colorant"yellow"), length=nr_colors)
colors = distinguishable_colors(nr_colors, 
    [RGB(1, 1, 1), RGB(0, 0, 0), RGB(0.8, 1.0, 0.8),RGB(0.6274, 1.0, 1.0), RGB(1.0, 1.0, 0.0)], 
    dropseed=true, lchoices=range(0, stop=80, length=15),)
idx_to_color = Dict(zip(idxs_active,colors))

pl_ts = plot( ylabel=L"L_i(t_L)", xlabel="", legend=false, xformatter=:none,
    tickfontsize=6, guidefontsize=18, left_margin=13Plots.mm)

#add vertical lines 
Δt_chart = 1.0
n_Δt = Int(Δt_chart / Δt)
vline!(pl_ts, t[tidxs[1]:n_Δt:tidxs[end]], ls=:dash, lc=:gray50, la=0.5)


#plot trajectory itself
pl_ts = plot!(pl_ts,t[tidxs], traj[tidxs, idxs_active], lc=colors', lw=3)
   


#plot inactive with gray separately
idxs_inactive = []
for i in 1:N
    i ∈ idxs_active && continue
    push!(idxs_inactive,i)
end

plot!(pl_ts,t[tidxs], traj[tidxs, idxs_inactive], lc=:gray50, lw=1,
    ylabel=L"L_i(t_L)", xlabel="", legend=false, xformatter=:none,
    tickfontsize=10, guidefontsize=18, left_margin=13Plots.mm)

#plot integrated variables
pl_int_act = plot(
  ylabel=L"\tilde{L}_i(t_L,\Delta t_L)")

#ylabel = L"\int_{Δt_i}^{Δt_{i+1}} L_i(t)Δt"
#add vertical lines 
vline!(pl_int_act, t_int[tint_idxs], ls=:dash, lc=:gray50, la=0.5)

plot!(pl_int_act,t_int[tint_idxs], int_act_matrix[tint_idxs, idxs_active], 
    lc=colors', mc=colors', lw=0.5, markershape=:circle, ms=6,
    markerstrokewiΔth=0.0, markerstrokecolor=colors')

#plot inactive with gray separately
plot!(t_int[tint_idxs], int_act_matrix[tint_idxs, idxs_inactive],
    lc=:gray50, mc=:gray50, lw=1,
    legend=false, tickfontsize=10, guidefontsize=15,xformatter=:none)

#plot charts and ranks
pl_topk = plot(legend=false,dpi=300,
    ylabel="ranks",xlabel=L"t_L",guidefontsize=10,
    left_margin=13Plots.mm,bottom_margin=4Plots.mm,
    yflip=true, yticks=[1:k;], ylims=(0, 5+0.5),
    tickfontsize=11)

vline!(pl_topk, t_int[tint_idxs], ls=:dash, lc=:gray50, la=0.5)

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


#add fig labels
x_annot = 7.73
annotate!(pl_ts, x_annot, 0.24, text(L"(a)", :gray10, :center, 18))
annotate!(pl_int_act, x_annot, 0.14, text(L"(b)", :gray10, :center, 18))
annotate!(pl_topk, x_annot, -0.2, text(L"(c)", :gray10, :center, 18))

plot!(pl_ts,ylims=(0.0,0.25),yticks=[0.0,0.05,0.1,0.15,0.2])
plot!(pl_int_act, ylims=(0.0, 0.15), yticks=[0.0, 0.05, 0.1])
plot!(pl_topk,xticks=[10:5:tspan[end];],guidefontsize=18)


pl = plot(pl_ts, pl_int_act, pl_topk, layout=(3, 1), size=(700, 400),dpi=300)

savefig(pl,"../../paper_figs/chart_explanation_example_shorter.png")