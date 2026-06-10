using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#--------------------------------plot N300 phase diagram------------------------------
N=300

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
data_dir_root = "../../data/phase_diag_N300_recalc_fixed/"
#data_dir_root = "../../data/phase_diag_N300_fp/"

dirs_pd300_fp = readdir(data_dir_root_fp)
dirs_pd300 = readdir(data_dir_root)
files_pd_levels = dirs_pd300[2:end-1] #[[2:17...,19:81...,18]]
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
    Λs = lyap_data_sorted[:, end] / 0.005
    Λs[Λs.<0] .= Λs[Λs.<0] * fp_scale #scale negative values

    ws_level = lyap_data_sorted[:, 3]
    push!(ws_vals, ws_level)
    phase_diag[i, 1:length(ws_level)] .= Λs


    #=
    #overwrite fixed point part 

    level_dir_fp = files_pd_levels_fp[i]
    data_dir_fp = data_dir_root_fp * level_dir_fp * "/"
    data_file_fp = readdir(data_dir_fp)[1]
    lyap_data_fp = readdlm(data_dir_fp * data_file_fp) 
    lyap_data_sorted_fp = lyap_data_fp[sortperm(lyap_data_fp[:,3]),:] 
    Λs_fp = lyap_data_sorted_fp[:, end] ./ 0.005 
    Λs_fp[Λs_fp.<0] .= Λs_fp[Λs_fp.<0] * fp_scale #scale negative values

    ws_fp_level = lyap_data_sorted_fp[:,3]
    phase_diag[i, 1:length(ws_fp_level)] .= Λs_fp
    =#

end

ws = ws_vals[1]

pl_pd300 = plot_phase_diagram(phase_diag, ws, g_rel_levels, N;
    xlims=(0.0, 0.95),
    ylims=(g_rel_levels[1], g_rel_levels[end]),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"1/(T_x(N-1))",clims=(-7,7),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar = false,
    legend=false,
    x_ticks = [0.0, 0.45, 0.9])


#add r curve 
r_curve(w, K, c) = K * c * (1 / w - 1) #g/(N-1) = 1/(ϵ*(N-1)), rp=rc
c0 = 2.4
K = 1.0 

#horizontal direction on phasediagram 
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c_article/" #"c/"
params_orig_horizontal = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled_horizontal = readdlm(param_dir * "params_w_g.txt")
params_changed_horizontal = readdlm(param_dir * "params_c.txt")
cs = params_changed_horizontal
ws_horizontal = params_rescaled_horizontal[:, 1]
gs_horizontal = params_rescaled_horizontal[:, 2]
g_horizontal = gs_horizontal[1]

param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_alpha_N.txt")
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = params_changed_curve[:, 1]
rs_rounded = round.(rs; digits=3)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2]

_, w_B2, w_B3 = ws_horizontal

w_B3,w_B2,w_B4 = ws_curve
ws_B3_B4 = [w_B3:0.001:w_B4;]
r_curve_vals = r_curve.(ws_B3_B4, K, c0)

#add r curve to plot 
plot!(ws_B3_B4, r_curve_vals, ls=:dot, lw=2, lc=:white, 
    ylims=(0, 10))#label=L" \propto (1/w - 1)")

g_B3, g_B2, g_B4 = gs_curve

#add points 
plot!([w_B2], [g_B2], ms=6, st=:scatter, mc=:lime,label=L"B_2",markershape=:circle)
plot!([w_B3], [g_B3], ms=6, st=:scatter, mc=:orange, label= L"B_3",markershape=:utriangle)
plot!([w_B4], [g_B4], ms=6, st=:scatter, mc=:purple, label = L"B_4",markershape=:rect)

#add two special levels 
hline!([g_B2],ls=:dash,lw=2,lc=:lime)
hline!([g_B4],ls=:dash,lw=2,lc=:purple)

#add annotations
annotate!(w_B3 + 0.08, g_B3 - 0.5, text(L"B_3", :orange, :center, annot_size))
annotate!(w_B2 + 0.08, g_B2 - 0.6, text(L"B_2", :lime, :center, annot_size))
annotate!(w_B4 + 0.04, g_B4 + 0.6, text(L"B_4", :purple, :center, annot_size))


#add colorbar 
cbar_range = [-7:0.1:7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-5,0.0,5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-6Plots.mm,tickfontsize=tick_fontsize_large-5)

#combine heatmap with the separate colorbar
l = @layout [a{0.95w} b]
pl_phase_diag = plot!(pl_pd300, pl_cbar, layout=l, tick_direction=:out)

#plot two levels of lyap curves separately
g_B4_approx = g_rel_levels[24]
g_B2_approx = g_rel_levels[53]

pl_curves = plot(ws,phase_diag[24,:],lc=:purple)
plot!(ws,phase_diag[53,:],lc=:lime)
plot!([w_B2], [g_B2],st=:scatter,ms=8,markershape=:circle)
plot!([w_B4], [g_B4],st=:scatter,ms=8,markershape=:rect)


#save figures and save processed data

savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_unfinished.png")
writedlm("../../data/phase_diag_N300/phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",phase_diag)
writedlm("../../data/phase_diag_N300/w_values_phase_diag_N300_Ttrmax_500_times_scaled.txt",ws)
writedlm("../../data/phase_diag_N300/g_rel_levels_phase_diag_N300_Ttrmax_500_times_scaled.txt",g_rel_levels)