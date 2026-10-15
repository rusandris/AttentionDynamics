using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
cd(@__DIR__)

include("../../plotting/plot_functions/lyapunov.jl")
include("../../plotting/paper_fig_plotters/defaults.jl")

r_curve(w, K, c) = K * c * (1 / w - 1)
N = 300
K0 = 1
c0 = 2.4
α0 = 0.005

src_dir = "../../data/phase_diag_N300_extended_highres/"
dirs_all = readdir(src_dir) 
g_rel_level_dirs = dirs_all[2:end] .* "/"
lyaps_data = []

for dir in g_rel_level_dirs
    f = readdir(src_dir * dir)[1]
    @show f
    λs = readdlm(src_dir * dir * f)
    push!(lyaps_data,λs)
end

phase_diag_vecs = []
for i in 1:length(lyaps_data)
    order_vals = sortperm(lyaps_data[i][:,4])
    λs = lyaps_data[i][:,end][order_vals]
    push!(phase_diag_vecs,λs)
end

g_rel_levels_ext = [data[1,1] for data in lyaps_data]
ws_ext = sort(lyaps_data[1][:,3])

phase_diag_ext = hcat(phase_diag_vecs...)
pd_scaled_ext = phase_diag_ext' ./ α0
pd_scaled_ext[pd_scaled_ext .< 0] .= pd_scaled_ext[pd_scaled_ext .< 0] .* 10

pl = plot_phase_diagram(pd_scaled_ext, ws_ext, g_rel_levels_ext, N;
    xlims=(0.0, 0.4),
    ylims=(g_rel_levels_ext[1], g_rel_levels_ext[end]),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"y",clims=(-10,10), #ylabel=L"1/(T_x(N-1))"
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar = false,
    legend=:topright,
    legend_fontsize=12)

plot!(pl,ws_ext,r_curve.(ws_ext,K0,c0),ls=:dot, lw=3, lc=:white,label="")

annotate!(pl,0.23, 14.5, text(L"r", :white, :center, 35),label="")
plot!(pl,[0.22,0.2],[13.5,15.0],arrow=:arrow,lc=:white,mc=:white,lw=4,label="")

#Books
rs_books = [8.0,12.0,15.0,17.0,18.0,20.0]
N_books = 300
K =1
resc_params_books = rescaled_params.(rs_books, rs_books, K, 2.4, 0.003, N_books)
ws_books = [resc_params_books[i][1] for i in 1:length(resc_params_books)]
norm_tscale_books = [resc_params_books[i][2] for i in 1:length(resc_params_books)] ./ (N_books - 1)


plot!(ws_books, norm_tscale_books, st=:scatter,markershape=:x,markerstrokewidth=6, ms=6, lc=:orange, mc=:orange, label="Google books")




#savefig(pl,"../../paper_figs/phase_diag_N300_extended_y10_20_T2000_Ttr2000.png")