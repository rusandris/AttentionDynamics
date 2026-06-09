using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

data_dir_root = "../../data/phase_diag_N300_fp_fixed/"
#data_dir_root = "../../data/phase_diag_N300_fp/"

dirs_pd300 = readdir(data_dir_root)
files_pd_levels = dirs_pd300[[2:17...,19:81...,18]]
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
#fill up phase_diag
for (i,level_dir) in enumerate(files_pd_levels)
    @show level_dir
    data_dir = data_dir_root * level_dir * "/"
    data_file = readdir(data_dir)[1]
    lyap_data = readdlm(data_dir * data_file) 
    lyap_data_sorted = lyap_data[sortperm(lyap_data[:,3]),:] 
    Λs = lyap_data_sorted[:,end] / 0.005
    ws_level = lyap_data_sorted[:,3]
    push!(ws_vals,ws_level)
    phase_diag[i ,1:length(ws_level)] .= Λs
end

c = cgrad(:RdBu, rev=true)
#c = cgrad([:blue,:white,:red], [0.05,0.1], categorical = true)
ws = ws_vals[1]
pl_pd300 = heatmap(ws,g_rel_levels,phase_diag .* 10, yflip=false,cmap=c,clims=(-8,8),
    dpi=300,xlabel=L"w",ylabel=L"1/(T_x(N-1))",guidefontsize=25,xticks=[0.1:0.1:1.0;],xlims=(0,1))

savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled.png")
writedlm("../../data/phase_diag_N300/phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",phase_diag)
writedlm("../../data/phase_diag_N300/w_values_phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",ws)
writedlm("../../data/phase_diag_N300/g_rel_levels_phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",g_rel_levels)