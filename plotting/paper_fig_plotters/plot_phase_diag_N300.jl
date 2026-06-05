using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

data_dir_root = "../../data/phase_diag_N300/"

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
ws = zeros(nr_ws)
#fill up phase_diag
for (i,level_dir) in enumerate(files_pd_levels)
    lyap_data = readdlm(data_dir_root * level_dir * "/lyaps_N300.txt") 
    lyap_data_sorted = lyap_data[sortperm(lyap_data[:,3]),:] 
    Λs = lyap_data_sorted[:,end] / 0.005
    ws .= lyap_data_sorted[:,3]
    length(Λs) < nr_ws && continue 

    @show length(Λs)
     phase_diag[i ,:] .= Λs
end
c = cgrad(:RdBu, [0.05,0.14], rev=true)
#c = cgrad([:blue,:white,:red], [0.05,0.1], categorical = true)
pl_pd300 = heatmap(ws,g_rel_levels,phase_diag, yflip=false,cmap=c,clims=extrema(phase_diag),
    dpi=300,xlabel=L"w",ylabel=L"1/(T_x(N-1))",guidefontsize=25,xticks=[0.1:0.1:0.95;])

savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled.png")