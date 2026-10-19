using Plots,LaTeXStrings
using DelimitedFiles
using AttentionDynamics
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

function read_phasediag_data(data_dir_root) 
    
    #read subdirs from root 
    subdirs = readdir(data_dir_root)
    y_values_dir = subdirs[findfirst(d -> occursin("g_rel_values",d),subdirs)]
    lyap_subdirs = filter(d -> d !== y_values_dir && occursin("lyapunovs",d),subdirs)
    sort!(lyap_subdirs;by=extract_num)
    @show lyap_subdirs

    #read y values
    files_y = readdir(data_dir_root * y_values_dir)
    y_levels_batches = [] 
    for file in files_y
        y_vals_batch = readdlm(data_dir_root * y_values_dir * "/" * file)
        push!(y_levels_batches,y_vals_batch) 
    end
    y_levels = vec(vcat(y_levels_batches...))
    nr_levels = length(y_levels)

    #init matrix for phase diagram
    probe_file = readdir(data_dir_root * lyap_subdirs[1])[1] 
    ws = sort(readdlm(data_dir_root * lyap_subdirs[1] * "/" * probe_file)[:,3])
    nr_ws = length(ws)
    @show nr_ws
    phase_diag = zeros(nr_levels,nr_ws)

    #fill up phase_diag
    for i in 1:length(y_levels)

        #skip unfinished values 
        i > nr_levels && continue     

        #fill lyaps
        level_dir = lyap_subdirs[i]
        data_dir = data_dir_root * level_dir * "/"
        files = readdir(data_dir)
        data_file = files[findfirst(f -> occursin("lyaps",f),files)]

        #find lyap file file 
        @show level_dir
        lyap_data = readdlm(data_dir * data_file)
        lyap_data_sorted = lyap_data[sortperm(lyap_data[:, 3]), :]    
    
        #scale lyap exp values
        Λs = lyap_data_sorted[:, end] / α0
        Λs[Λs.<0] .= Λs[Λs.<0] * fp_scale #scale negative values
        phase_diag[i, 1:length(ws)] .= Λs

    end
    return ws,y_levels,phase_diag
end 

function extract_num(x::String)
    number_idx_start = findfirst(isdigit, x)
    number_idx_stop = findnext(c-> !(isdigit(c)) && c != '.', x,number_idx_start)-1
    num = SubString(x, number_idx_start,number_idx_stop)
    return parse(Float64, num)
end

r_curve(w, K, c) = K * c * (1 / w - 1) #g/(N-1) = 1/(ϵ*(N-1)), rp=rc


#--------------------------------plot N300 phase diagram------------------------------
N=300
α0 = 0.005
K = 1.0
c0 = 2.4

data_dir_root_fp = "../../data/phase_diag_N300_fp_fixed/"
dirs_pd300_fp = readdir(data_dir_root_fp)

#phasediag data: y in [0,10]
data_dir_root = "../../data/phase_diag_N300_recalc_fixed_T500_Ttr500/"
dirs_pd300_all = readdir(data_dir_root)


#phasediag data: y in [10,20]
data_extension_dir_root = "../../data/phase_diag_N300_extended_sameres_matched/"
#data_dir_root = "../../data/phase_diag_N300_recalc_fixed_T2000_Ttr2000//"
#data_dir_root = "../../data/phase_diag_N300_fp/"

#data constants 
fp_scale = 10
#ws_fp_idx_offset = 9

#----------------------------------------------Fig4-----------------------------------------------

ws,ys,pd = read_phasediag_data(data_dir_root)


pl_pd300 = plot_phase_diagram(pd, ws, ys, N;
    xlims=(0.0, 0.95),
    ylims=(ys[1], ys[end]),
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
    legend_fontsize=12,
    x_ticks = [0,0.45, 0.9])

#make a copy before adding stuff 
pl_pd300_zoom_ext = deepcopy(pl_pd300)

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

#--------------------------------add markers and curves to phase diagram----------------------------------

colors = [:orange,:cyan,:slateblue1,:gray70,:royalblue1,:lightseagreen,:lime]
ms_data = 5
plot!(pl_pd300,ws_B2_B4, r_curve_vals, ls=:dot, lw=3, lc=:white, ylims=(0,10),label="") #label=L" \propto (1/w - 1)")


#add special points with special markers to phasediag----------------------
ms = 8
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

#savefig(pl_pd300,"../../paper_figs/phase_diag_N300_Ttrmax_500_scaled_times.png")
savedir = "../../data/phase_diag_N300_recalc_fixed/process_results/"
mkpath(savedir)
#=
writedlm(savedir * "phase_diag_N300_Ttrmax_500_times_scaled_matrix_form.txt",phase_diag)
writedlm(savedir * "w_values_phase_diag_N300_Ttrmax_500_times_scaled.txt", ws)
writedlm(savedir * "g_rel_levels_phase_diag_N300_Ttrmax_500_times_scaled.txt", g_rel_levels)
=#

#calculate (y,w) rescaled parameters for every parameter for every dataset -> df::DataFrame 
include("generate_parameters_datasets_table.jl")

#plot every data set

msw = 5
mshape = :x

#Books dataset parameters 
df_books = df[occursin.("Books",df.Dataset),:] 
ws_books = df_books.w
y_books = df_books.y 

#Twitter dataset parameters 
df_twitter = df[occursin.("Twitter",df.Dataset),:] 
ws_twitter = df_twitter.w
y_twitter = df_twitter.y 

#Publications dataset parameters 
df_pub = df[occursin.("Publications",df.Dataset),:] 
ws_pub = df_pub.w
y_pub = df_pub.y 

#movies dataset parameters 
df_mov = df[occursin.("Movies",df.Dataset),:] 
ws_mov = df_mov.w
y_mov = df_mov.y 

#Twitter dataset parameters 
df_reddit = df[occursin.("Reddit",df.Dataset),:] 
ws_reddit = df_reddit.w
y_reddit = df_reddit.y 

#Publications dataset parameters 
df_trends = df[occursin.("Publications",df.Dataset),:] 
ws_trends = df_trends.w
y_trends = df_trends.y 


#leave out wiki
pl_every_data = plot([1e-4:1e-4:0.4;], 1 ./ [1e-4:1e-4:0.4;], ls=:dash, lw=3, lc=:gray10, label="",yscale=:log10)
plot!(pl_every_data,[1e-4:1e-4:0.4;], r_curve.([1e-4:1e-4:0.4;], K, c0), ls=:dot, lw=3, lc=:gray70, mc=:gray10, label="")
plot!(pl_every_data,ws_books, y_books, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[1], mc=colors[1], label="Google books")
plot!(pl_every_data,ws_mov, y_mov,st=:scatter ,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[2], mc=colors[2], label="Movies")
plot!(pl_every_data,ws_trends, y_trends, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[3], mc=colors[3], label="Google Trends")
plot!(pl_every_data,ws_reddit, y_reddit, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[4], mc=colors[4], label="Reddit")
plot!(pl_every_data,ws_pub, y_pub,st=:scatter ,markershape=mshape,markerstrokewidth=msw, ms=ms_data+3, lc=colors[5], mc=colors[5], label="Publications")
#plot!(ws_wiki, norm_tscale_wiki, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[6], mc=colors[6], label="Wikipedia")
plot!(pl_every_data,ws_twitter, y_twitter, st=:scatter,markershape=mshape,markerstrokewidth=msw, ms=ms_data, lc=colors[7], mc=colors[7], label="Twitter")

plot!(pl_every_data,legendfontsize=legend_fontsize-5, guidefontsize=25, tickfontsize=15, 
    xlabel=L"w", ylabel=L"y",dpi=300,xlims=(-0.01,0.4),size=(600,400),
    xticks=[0,0.1,0.2,0.3,0.4],yticks=[0,1e1,1e2,1e3],ylims=(1.0,1e3),
    top_margin=5Plots.mm)

#savefig(pl_every_data,"../../paper_figs/phase_diag_all_data_curves_SupFig10_log_yscale.png")


#----------------------------------------------Fig3-----------------------------------------------

#---------------------------combine y < 10 phasediag and y > 10 phasediag------------------------------

ws_ext,ys_ext,pd_ext = read_phasediag_data(data_extension_dir_root)

#crop y < 10 phase diag
phase_diag_zoom_ext = vcat(pd[:,1:length(ws_ext)],pd_ext[2:end,:])
ys_zoom_ext = vcat(ys,ys_ext[2:end])


pl_pd300_zoom_ext = plot_phase_diagram(phase_diag_zoom_ext, ws_ext, ys_zoom_ext, N;
    xlims=(0.0, 0.4),
    ylims=(ys_zoom_ext[1], ys_zoom_ext[end]),
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
    legend_fontsize=12,
    xticks=[0,0.1,0.2,0.3,0.4])


#add r curve to plot 
ws_curve = [0.1:0.001:0.4;]
r_curve_vals = r_curve.(ws_curve, K, c0)
plot!(ws_curve, r_curve_vals, ls=:dot, lw=3, lc=:white, ylims=(0,20),label="") #label=L" \propto (1/w - 1)")


#Books, Publications, Twitter
plot!(pl_pd300_zoom_ext,ws_curve, r_curve_vals, ls=:dot, lw=3, lc=:white, ylims=(0,20),label="") #label=L" \propto (1/w - 1)")
plot!(pl_pd300_zoom_ext, ws_books, y_books, st=:scatter,ms=ms_data, mc=colors[1],markershape=:x,
    markerstrokewidth=msw, ylims=(0,20), label="Google books",legend=:bottomleft,legendfontsize = legend_fontsize-5) #label=L" \propto (1/w - 1)")
plot!(pl_pd300_zoom_ext, ws_pub, y_pub, st=:scatter, ms=ms_data+3, mc=colors[5],markershape=:x,markerstrokewidth=msw, ylims=(0, 20), label="Publications") #label=L" \propto (1/w - 1)")
plot!(pl_pd300_zoom_ext, ws_twitter, y_twitter, st=:scatter, ms=ms_data,markershape=:x,markerstrokewidth=msw, mc=colors[7], ylims=(0, 20), label="Twitter") #label=L" \propto (1/w - 1)")

#add twitter year annotations
annot_offset_x = 0.02
annot_offset_y = 0.1
annotate!(pl_pd300_zoom_ext,ws_twitter[1]+annot_offset_x, y_twitter[1]+annot_offset_y, text(L"2013", :lime, :left, annot_size_letters-2))
annotate!(pl_pd300_zoom_ext,ws_twitter[2]+annot_offset_x, y_twitter[2]+annot_offset_y, text(L"2014", :lime, :left, annot_size_letters-2))
annotate!(pl_pd300_zoom_ext,ws_twitter[3]+annot_offset_x, y_twitter[3]+annot_offset_y, text(L"2015", :lime, :left, annot_size_letters-2))
annotate!(pl_pd300_zoom_ext,ws_twitter[4]+annot_offset_x, y_twitter[4]+annot_offset_y, text(L"2016", :lime, :left, annot_size_letters-2))

#plot!([w_B1], [g_B2], ms=8, st=:scatter, mc=:cyan,label=L"B_1",markershape=:diamond)
plot!(pl_pd300_zoom_ext,[w_B2], [g_B2], ms=8, st=:scatter, mc=:lime,label="",markershape=:dtriangle)
plot!(pl_pd300_zoom_ext,[w_B3], [g_B3], ms=8, st=:scatter, mc=:orange, label= "",markershape=:utriangle)
plot!(pl_pd300_zoom_ext,[w_B4], [g_B4], ms=8, st=:scatter, mc=:purple, label = "",markershape=:rect)

annotate!(pl_pd300_zoom_ext,w_B3 - 0.01, g_B3 - 0.3, text(L"B_3", :orange, :right, annot_size_letters))
annotate!(pl_pd300_zoom_ext,w_B2 - 0.01, g_B2 - 0.1, text(L"B_2", :lime, :right, annot_size_letters))

#combine phase diag with data and log scale all data fig 

pl_phase_diag_simple = plot(pl_pd300_zoom_ext, pl_cbar, layout=l_cbar, tick_direction=:out)
#add cbar annot
annotate!(pl_phase_diag_simple, 0.418, 10.5, text(L"\Lambda", :gray10, :center, annot_size))

pl_fig3_v1 = plot(pl_phase_diag_simple,pl_every_data,layout=(1,2),size=(1000,400),
    bottom_margin=6Plots.mm,top_margin=10Plots.mm) 

#=
#add fig level annotation 
annotate!(pl_fig3_v1[1],-0.08, 10.8, text(L"(a)", :gray10, :center, annot_size))
annotate!(pl_fig3_v1[3], -0.1, 17.0*1e2, text(L"(b)", :gray10, :center, annot_size))
plot!(pl_fig3_v1[2],right_margin=6Plots.mm) #adjust right margin
plot!(pl_fig3_v1[3],right_margin=6Plots.mm) #adjust right margin
=#

#savefig(pl_fig3_v1,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_all_data_curves_SupFig10_log_yscale.png")


