using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
data_dir_root = "../../data/phase_diag_N300/"
#data_dir_root = "../../data/phase_diag_N300_fp/"

dirs_pd300_fp = readdir(data_dir_root_fp)
dirs_pd300 = readdir(data_dir_root)
files_pd_levels = dirs_pd300[[2:17...,19:81...,18]]
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
phase_diag = zeros(nr_levels,nr_ws)
ws_vals = Vector{Float64}[]
ws_fp_idx_offset = 9
#fill up phase_diag
for i in 1:length(g_rel_levels)

    #fill lyaps  
    level_dir = files_pd_levels[i]
    data_dir = data_dir_root * level_dir * "/"
    data_file = readdir(data_dir)[1]
    @show level_dir
    lyap_data = readdlm(data_dir * data_file)
    lyap_data_sorted = lyap_data[sortperm(lyap_data[:, 3]), :]
    Λs = lyap_data_sorted[:, end] / 0.005
    ws_level = lyap_data_sorted[:, 3]
    push!(ws_vals, ws_level)
    phase_diag[i, 1:length(ws_level)] .= Λs


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

end

c = cgrad(:RdBu, rev=true)
#c = cgrad([:blue,:white,:red], [0.05,0.1], categorical = true)
ws = ws_vals[1]
pl_pd300 = heatmap(ws,g_rel_levels,phase_diag, yflip=false,cmap=c,clims=(-7,7),
    dpi=300,xlabel=L"w",ylabel=L"1/(T_x(N-1))",guidefontsize=25,xticks=[0.1:0.1:1.0;],xlims=(0,1))

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

c = g_B2 / K / (1 / w_B2 - 1)

w_B3,w_B2,w_B4 = ws_curve
ws_B3_B4 = [w_B3:0.001:w_B4;]
r_curve_vals = r_curve.(ws_B3_B4, K, c0)

#add r curve to plot 
plot!(ws_B3_B4, r_curve_vals, ls=:dash, lw=1, lc=:gray70, 
    ylims=(0, 10),label=L" \propto (1/w - 1)")

g_B2 = g_horizontal
g_B3 = gs_alpha[1]
g_B4 = gs_alpha[end]

#add points 
plot!([w_B2], [g_B2], ms=6, st=:scatter, mc=:lime,label=L"B_2")
plot!([w_B3], [g_B3], ms=6, st=:scatter, mc=:orange, label= L"B_3")
plot!([w_B4], [g_B4], ms=6, st=:scatter, mc=:purple, label = L"B_4")
plot!(legendfontsize=12,legend=:topright)

savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_combined_artifical.png")
writedlm("../../data/phase_diag_N300/phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",phase_diag)
writedlm("../../data/phase_diag_N300/w_values_phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",ws)
writedlm("../../data/phase_diag_N300/g_rel_levels_phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",g_rel_levels)