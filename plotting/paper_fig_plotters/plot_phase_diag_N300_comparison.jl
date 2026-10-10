using Plots,LaTeXStrings
using DelimitedFiles
using AttentionDynamics
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#--------------------------------plot N300 phase diagram------------------------------

#-----------------------------------T=Ttr=2000------------------------
N=300
α0 = 0.005

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
data_dir_root = "../../data/phase_diag_N300_recalc_fixed_T2000_Ttr2000//"
#data_dir_root = "../../data/phase_diag_N300_recalc_fixed_T500_Ttr500//"
#data_dir_root = "../../data/phase_diag_N300_fp/"

dirs_pd300_fp = readdir(data_dir_root_fp)
dirs_pd300 = readdir(data_dir_root)
files_pd_levels = dirs_pd300[[2:17...,19:81...,18]] #[2:end] #
files_pd_levels_fp = dirs_pd300_fp[[2:17..., 19:81..., 18]]

fp_scale = 10
nr_levels = length(files_pd_levels)
nr_ws = 80

#read g_rel_values
files_grel = readdir(data_dir_root * dirs_pd300[1])
g_rel_levels_batches = [] 
for file in files_grel
    g_vals_batch = readdlm(data_dir_root * dirs_pd300[1] * "/" * file)
    push!(g_rel_levels_batches,g_vals_batch) 
end
g_rel_levels = vec(vcat(g_rel_levels_batches...))

#init matrix for phase diagram 
phase_diag = zeros(nr_ws,nr_ws)
ws_vals = Vector{Float64}[]
ws_fp_idx_offset = 9

#fill up phase_diag
for i in 1:length(g_rel_levels)

    #skip unfinished values 
    i > nr_levels && continue     

    #fill lyaps  
    level_dir = files_pd_levels[i]
    data_dir = data_dir_root * level_dir * "/"
    data_file = readdir(data_dir)[1]
    @show level_dir
    lyap_data = readdlm(data_dir * data_file)
    lyap_data_sorted = lyap_data[sortperm(lyap_data[:, 3]), :]
    Λs = lyap_data_sorted[:, end] / α0
    Λs[Λs.<0] .= Λs[Λs.<0] * fp_scale #scale negative values

    ws_level = lyap_data_sorted[:, 3]
    push!(ws_vals, ws_level)
    phase_diag[i, 1:length(ws_level)] .= Λs

end

ws = ws_vals[1]

pl_pd300_T2000 = plot_phase_diagram(phase_diag, ws, g_rel_levels, N;
    xlims=(0.0, 0.95),
    ylims=(g_rel_levels[1], g_rel_levels[end]),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"y",clims=(-7,7), #ylabel=L"1/(T_x(N-1))"
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar = false,
    legend=:topright,
    legend_fontsize=12,
    x_ticks = [0,0.45, 0.9])

plot!(pl_pd300_T2000,right_margin=10Plots.mm,title=L"T=T_{tr}=2000")

#add colorbar 
cbar_range = [-7:0.1:7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-5,0.0,5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-8Plots.mm,tickfontsize=tick_fontsize_large)

#combine heatmap with the separate colorbar
l_cbar = @layout [a{0.95w} b]
pl_phase_diag_T2000 = plot!(pl_pd300_T2000, pl_cbar, layout=l_cbar, tick_direction=:out)

#add cbar annot
annotate!(pl_phase_diag_T2000,1.0, 10.5, text(L"\Lambda", :gray10, :center, annot_size))


#-----------------------------------T=Ttr=500------------------------
N=300
α0 = 0.005

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
data_dir_root = "../../data/phase_diag_N300_recalc_fixed_T500_Ttr500//"
#data_dir_root = "../../data/phase_diag_N300_fp/"

dirs_pd300_fp = readdir(data_dir_root_fp)
dirs_pd300 = readdir(data_dir_root)
files_pd_levels = dirs_pd300[[2:17...,19:81...,18]] #[2:end] #
files_pd_levels_fp = dirs_pd300_fp[[2:17..., 19:81..., 18]]

fp_scale = 10
nr_levels = length(files_pd_levels)
nr_ws = 80

#read g_rel_values
files_grel = readdir(data_dir_root * dirs_pd300[1])
g_rel_levels_batches = [] 
for file in files_grel
    g_vals_batch = readdlm(data_dir_root * dirs_pd300[1] * "/" * file)
    push!(g_rel_levels_batches,g_vals_batch) 
end
g_rel_levels = vec(vcat(g_rel_levels_batches...))

#init matrix for phase diagram 
phase_diag = zeros(nr_ws,nr_ws)
ws_vals = Vector{Float64}[]
ws_fp_idx_offset = 9

#fill up phase_diag
for i in 1:length(g_rel_levels)

    #skip unfinished values 
    i > nr_levels && continue     

    #fill lyaps  
    level_dir = files_pd_levels[i]
    data_dir = data_dir_root * level_dir * "/"
    data_file = readdir(data_dir)[1]
    @show level_dir
    lyap_data = readdlm(data_dir * data_file)
    lyap_data_sorted = lyap_data[sortperm(lyap_data[:, 3]), :]
    Λs = lyap_data_sorted[:, end] / α0
    Λs[Λs.<0] .= Λs[Λs.<0] * fp_scale #scale negative values

    ws_level = lyap_data_sorted[:, 3]
    push!(ws_vals, ws_level)
    phase_diag[i, 1:length(ws_level)] .= Λs

end

ws = ws_vals[1]

pl_pd300_T500 = plot_phase_diagram(phase_diag, ws, g_rel_levels, N;
    xlims=(0.0, 0.95),
    ylims=(g_rel_levels[1], g_rel_levels[end]),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"y",clims=(-7,7), #ylabel=L"1/(T_x(N-1))"
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar = false,
    legend=:topright,
    legend_fontsize=12,
    x_ticks = [0,0.45, 0.9])


plot!(pl_pd300_T500,left_margin=10Plots.mm,title=L"T=T_{tr}=500")


#add colorbar 
cbar_range = [-7:0.1:7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-5,0.0,5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-8Plots.mm,tickfontsize=tick_fontsize_large)


#combine heatmap with the separate colorbar
l_cbar = @layout [a{0.95w} b]
pl_phase_diag_T500 = plot!(pl_pd300_T500, pl_cbar, layout=l_cbar, tick_direction=:out)
#add cbar annot
annotate!(pl_phase_diag_T500,1.0, 10.5, text(L"\Lambda", :gray10, :center, annot_size))



pl = plot(pl_phase_diag_T500,pl_phase_diag_T2000,
    layout=(1,2),size=(1500,500),
    bottom_margin=10Plots.mm,
    top_margin=5Plots.mm)

savefig(pl,"../../paper_figs/phase_diag_N300_lyaps_log_comparison_T500_Ttr500_T2000_Ttr2000.png")