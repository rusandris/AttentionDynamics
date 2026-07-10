using Plots,LaTeXStrings
using DelimitedFiles
using AttentionDynamics
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#--------------------------------plot N300 phase diagram------------------------------
N=300
α0 = 0.005

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
data_dir_root = "../../data/phase_diag_N300_recalc_fixed/"
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

#make a copy before adding stuff 
pl_pd300_simple = deepcopy(pl_pd300)

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

w_B1, w_B2, w_B3 = ws_horizontal
g_B3, g_B2, g_B4 = gs_curve

w_B3,w_B2,w_B4 = ws_curve
ws_B2_B4 = [w_B2-0.1:0.001:w_B4+0.1;]
r_curve_vals = r_curve.(ws_B2_B4, K, c0)

#add r curve to plot 
plot!(ws_B2_B4, r_curve_vals, ls=:dot, lw=3, lc=:white, ylims=(0,10),label="") #label=L" \propto (1/w - 1)")
#writedlm("r_curve_values.txt",hcat(ws_B3_B4, r_curve_vals))

#TODO: try to add some data fit curves from Suppl. Fig10 from Lorenz-Spreen
# get parameters from captions, translate every parameter combination to (Tx,w)
# N=100 where does that fall on the N=300 phase diag?
colors = [:orange,:cyan,:slateblue1,:gray70,:royalblue1,:lightseagreen,:lime]
ms_data = 5
#Google Books
rs_books = [8.0,12.0,15.0,17.0,18.0,20.0]
N_books = 300 
resc_params_books = rescaled_params.(rs_books, rs_books, K, 2.4, 0.003, N_books)
ws_books = [resc_params_books[i][1] for i in 1:length(resc_params_books)]
norm_tscale_books = [resc_params_books[i][2] for i in 1:length(resc_params_books)] ./ (N_books - 1)

#Movies
rs_mov = [18.0, 20.0, 22.0, 25.0, 30.0, 33.0, 35.0, 37.0]
N_mov = 100
resc_params_mov = rescaled_params.(rs_mov, rs_mov, K, 2.4, 0.001, N_mov)
ws_mov = [resc_params_mov[i][1] for i in 1:length(resc_params_mov)]
norm_tscale_mov = [resc_params_mov[i][2] for i in 1:length(resc_params_mov)] ./ (N_mov - 1)

#Google Trends
rs_trends = [15.0, 17.0, 19.0, 21.0, 22.0, 26.0, 27.0, 30.0]
N_trends = 100
resc_params_trends = rescaled_params.(rs_trends, rs_trends, K, 2.0, 0.0005, N_trends)
ws_trends = [resc_params_trends[i][1] for i in 1:length(resc_params_trends)]
norm_tscale_trends = [resc_params_trends[i][2] for i in 1:length(resc_params_trends)] ./ (N_trends - 1)

#Reddit 
rs_reddit = [25.0, 29.0, 30.0, 33.0, 34.0, 35.0]
N_reddit = 100
resc_params_reddit = rescaled_params.(rs_reddit, rs_reddit, K, 5.4, 0.001, N_reddit)
ws_reddit = [resc_params_reddit[i][1] for i in 1:length(resc_params_reddit)]
norm_tscale_reddit = [resc_params_reddit[i][2] for i in 1:length(resc_params_reddit)] ./ (N_reddit - 1)

#publications 
rs_pub = [6.0, 7.0]
N_pub = 300
resc_params_pub = rescaled_params.(rs_pub, rs_pub, K, 2.4, 0.003, N_pub)
ws_pub = [resc_params_pub[i][1] for i in 1:length(resc_params_pub)]
norm_tscale_pub = [resc_params_pub[i][2] for i in 1:length(resc_params_pub)] ./ (N_pub - 1)

#wiki 
rs_wiki = [45.0, 48.0, 50.0]
N_wiki = 100
resc_params_wiki = rescaled_params.(rs_wiki, rs_wiki, K, 2.4, 0.0001, N_wiki)
ws_wiki = [resc_params_wiki[i][1] for i in 1:length(resc_params_wiki)]
norm_tscale_wiki = [resc_params_wiki[i][2] for i in 1:length(resc_params_wiki)] ./ (N_wiki - 1)

#Twitter
rs_twitter = [9.0,10.0,11.0,12.0]
N_twitter = 300 
resc_params_twitter = rescaled_params.(rs_twitter, rs_twitter, K, 2.4, 0.005, N_twitter)
ws_twitter = [resc_params_twitter[i][1] for i in 1:length(resc_params_twitter)]
norm_tscale_twitter = [resc_params_twitter[i][2] for i in 1:length(resc_params_twitter)] ./ (N_twitter - 1)


#edit simple phase diag with data points 
msw = 3
mshape = :x
plot!(pl_pd300_simple,ws_B2_B4, r_curve_vals, ls=:dot, lw=3, lc=:white, ylims=(0,10),label="") #label=L" \propto (1/w - 1)")
plot!(pl_pd300_simple, ws_books, norm_tscale_books, st=:scatter,ms=ms_data, mc=colors[1],markershape=:x,
    markerstrokewidth=msw, ylims=(0,10), label="Google books",legend=:bottomleft,legendfontsize = legend_fontsize-5) #label=L" \propto (1/w - 1)")
plot!(pl_pd300_simple, ws_pub, norm_tscale_pub, st=:scatter, ms=ms_data, mc=colors[5],markershape=:x,markerstrokewidth=msw, ylims=(0, 10), label="Publications") #label=L" \propto (1/w - 1)")
plot!(pl_pd300_simple, ws_twitter, norm_tscale_twitter, st=:scatter, ms=ms_data,markershape=:x,markerstrokewidth=msw, mc=colors[7], ylims=(0, 10), label="Twitter") #label=L" \propto (1/w - 1)")

annot_offset_x = 0.02
annot_offset_y = 0.1
annotate!(pl_pd300_simple,ws_twitter[1]+annot_offset_x, norm_tscale_twitter[1]+annot_offset_y, text(L"2013", :lime, :left, annot_size_letters-4))
annotate!(pl_pd300_simple,ws_twitter[2]+annot_offset_x, norm_tscale_twitter[2]+annot_offset_y, text(L"2014", :lime, :left, annot_size_letters-4))
annotate!(pl_pd300_simple,ws_twitter[3]+annot_offset_x, norm_tscale_twitter[3]+annot_offset_y, text(L"2015", :lime, :left, annot_size_letters-4))
annotate!(pl_pd300_simple,ws_twitter[4]+annot_offset_x, norm_tscale_twitter[4]+annot_offset_y, text(L"2016", :lime, :left, annot_size_letters-4))

#plot!([w_B1], [g_B2], ms=8, st=:scatter, mc=:cyan,label=L"B_1",markershape=:diamond)
plot!(pl_pd300_simple,[w_B2], [g_B2], ms=8, st=:scatter, mc=:lime,label="",markershape=:dtriangle)
plot!(pl_pd300_simple,[w_B3], [g_B3], ms=8, st=:scatter, mc=:orange, label= "",markershape=:utriangle)
plot!(pl_pd300_simple,[w_B4], [g_B4], ms=8, st=:scatter, mc=:purple, label = "",markershape=:rect)

annotate!(pl_pd300_simple,w_B3 - 0.01, g_B3 - 0.3, text(L"B_3", :orange, :right, annot_size_letters))
annotate!(pl_pd300_simple,w_B2 - 0.01, g_B2 - 0.1, text(L"B_2", :lime, :right, annot_size_letters))


#fall outside of phase diagram
#plot!(pl_pd300, ws_wiki, norm_tscale_wiki, st=:scatter, ms=ms_data, mc=colors[6], ylims=(0, 10), label="Wikipedia") #label=L" \propto (1/w - 1)")
#plot!(pl_pd300, ws_mov, norm_tscale_mov, st=:scatter,ms=ms_data, mc=colors[2], ylims=(0,10), label="Movies") #label=L" \propto (1/w - 1)")
#plot!(pl_pd300, ws_trends, norm_tscale_trends, st=:scatter, ms=ms_data, mc=colors[3], ylims=(0,10), label="Google trends") #label=L" \propto (1/w - 1)")
#plot!(pl_pd300, ws_reddit, norm_tscale_reddit, st=:scatter, ms=ms_data, mc=colors[4], ylims=(0,10), label="Reddit") #label=L" \propto (1/w - 1)")



#savefig(pl_pd300, "../../paper_figs/phase_diag_lyaps_N300_data_curves_SupFig10.png")

ms = 8

#add points 
#plot!([w_B1], [g_B2], ms=8, st=:scatter, mc=:cyan,label=L"B_1",markershape=:diamond)
plot!(pl_pd300,[w_B2], [g_B2], ms=ms, st=:scatter, mc=:lime,label="",markershape=:dtriangle)
plot!(pl_pd300,[w_B3], [g_B3], ms=ms, st=:scatter, mc=:orange, label= "",markershape=:utriangle)
plot!(pl_pd300,[w_B4], [g_B4], ms=ms, st=:scatter, mc=:purple, label = "",markershape=:rect)

#add two/three special levels 
hline!(pl_pd300,[g_B2],ls=:dash,lw=2,lc=:lime,label="")
hline!(pl_pd300,[g_B4],ls=:dash,lw=2,lc=:purple,label="")
hline!(pl_pd300,[g_B3],ls=:dash,lw=2,lc=:orange,label="")


#add annotations
annot_size_letters = annot_size - 3
#annotate!(w_B1 + 0.03, g_B2 + 0.6, text(L"B_1", :cyan, :center, annot_size_letters))
annotate!(pl_pd300,w_B3 + 0.05, g_B3 - 0.6, text(L"B_3", :orange, :center, annot_size_letters))
annotate!(pl_pd300,w_B2 + 0.05, g_B2 + 0.6, text(L"B_2", :lime, :center, annot_size_letters))
annotate!(pl_pd300,w_B4 + 0.04, g_B4 + 0.6, text(L"B_4", :purple, :center, annot_size_letters))

#add colorbar 
cbar_range = [-7:0.1:7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-5,0.0,5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-8Plots.mm,tickfontsize=tick_fontsize_large)

#combine heatmap with the separate colorbar
l_cbar = @layout [a{0.95w} b]
pl_phase_diag = plot!(pl_pd300, pl_cbar, layout=l_cbar, tick_direction=:out)

#add cbar annot
annotate!(pl_phase_diag,1.0, 10.5, text(L"\Lambda", :gray10, :center, annot_size))


#save figures and save processed data

savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_scaled_times.png")
savedir = "../../data/phase_diag_N300_recalc_fixed/process_results/"
mkpath(savedir)

writedlm(savedir * "phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",phase_diag)
writedlm(savedir * "w_values_phase_diag_N300_Ttrmax_500_times_scaled.txt", ws)
writedlm(savedir * "g_rel_levels_phase_diag_N300_Ttrmax_500_times_scaled.txt", g_rel_levels)

#plot every data set 
pl_every_data = plot([1e-4:1e-4:0.4;], 1 ./ [1e-4:1e-4:0.4;], ls=:dash, lw=3, lc=:gray10, label="",yscale=:log10)
plot!([1e-4:1e-4:0.4;], r_curve.([1e-4:1e-4:0.4;], K, c0), ls=:dot, lw=3, lc=:gray70, mc=:gray10, label="")
plot!(ws_books, norm_tscale_books, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[1], mc=colors[1], label="Google books")
plot!(ws_mov, norm_tscale_mov,st=:scatter ,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[2], mc=colors[2], label="Movies")
plot!(ws_trends, norm_tscale_trends, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[3], mc=colors[3], label="Google Trends")
plot!(ws_reddit, norm_tscale_reddit, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[4], mc=colors[4], label="Reddit")
plot!(ws_pub, norm_tscale_pub,st=:scatter ,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[5], mc=colors[5], label="Publications")
plot!(ws_wiki, norm_tscale_wiki, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[6], mc=colors[6], label="Wikipedia")
plot!(ws_twitter, norm_tscale_twitter, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[7], mc=colors[7], label="Twitter")

plot!(pl_every_data,legendfontsize=legend_fontsize-5, guidefontsize=25, tickfontsize=15, 
    xlabel=L"w", ylabel=L"y",dpi=300,xlims=(-0.01,0.4),size=(600,400),xticks=[0,0.1,0.2,0.3,0.4],yticks=[0,1e1,1e3])

savefig(pl_every_data,"../../paper_figs/phase_diag_all_data_curves_SupFig10_log_yscale.png")

#combine phase diag with data and log scale all data fig 
plot!(pl_pd300_simple, xticks=[0,0.1,0.2,0.3,0.4], xlims=(0.0,0.4))

pl_phase_diag_simple = plot(pl_pd300_simple, pl_cbar, layout=l_cbar, tick_direction=:out)
#add cbar annot
annotate!(pl_phase_diag_simple, 0.418, 10.5, text(L"\Lambda", :gray10, :center, annot_size))


pl_fig3_v1 = plot(pl_phase_diag_simple,pl_every_data,layout=(1,2),size=(1000,400),
    bottom_margin=6Plots.mm,top_margin=10Plots.mm) 

#add fig level annotation 
annotate!(pl_fig3_v1[1],-0.08, 10.8, text(L"(a)", :gray10, :center, annot_size))
annotate!(pl_fig3_v1[3], -0.08, 6.8*1e4, text(L"(b)", :gray10, :center, annot_size))
plot!(pl_fig3_v1[2],right_margin=6Plots.mm) #adjust right margin
plot!(pl_fig3_v1[3],right_margin=6Plots.mm) #adjust right margin


savefig(pl_fig3_v1,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_all_data_curves_SupFig10_log_yscale.png")


