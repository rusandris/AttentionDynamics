using Plots, LaTeXStrings
using AttentionDynamics
using DelimitedFiles
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#plot settings 
ms = 6
msw = 1.0
logocolors = Colors.JULIA_LOGO_COLORS
#-------------------------plot phase diagram---------------------

#--------------read lyaps--------------------
N = 10
figs_dir = "../../paper_figs/"
#data_dir = "../../data/phasediag_data_N10_extended_d0_1.0e-9_ut1e-6/"
#data_dir = "../../data/phasediag_data_N10_extended3_80x80/"
data_dir = "../../data/phasediag_data_N10_extended5_80x80/"
#lyap_vals = reshape(readdlm(data_dir*"lyap_values_Ttr10000.txt"),(80,80))
lyap_vals = reshape(readdlm(data_dir * "lyap_values_Ttr10000.txt"), (80, 80))
params = readdlm(data_dir*"phasediag_params.txt")
ws_phasediag = params[:,1]
gs_phasediag = params[:, 2] ./ (N-1)
N=10
chaos_curve(w) = w^(-1.5)


pl_lyaps = plot_phase_diagram(lyap_vals, ws_phasediag, gs_phasediag, N;
    xlims=(0.0, 1.0),
    ylims=(gs_phasediag[1], gs_phasediag[end]),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"1/(T_x(N-1))",clims=(-1,1),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",legend=false,
    colorbar = false,
    x_ticks = [0.0, 0.5, 1.0])


#plot!(pl_lyaps, ws_phasediag, chaos_curve.(ws_phasediag),ls=:dash,lw=2,lc=:purple)

#---------------------------fill in with analytical eigenvalues-------------------------------
λa = -1
λ_b(g,w,N) = g*(1-1/(N-1)*1/(1/w-1)) 
λ_d_realpart = -0.5

lyap_vals_filled = deepcopy(lyap_vals)
for (j,g) in enumerate(gs_phasediag) 
    for (i, w) in enumerate(ws_phasediag)
        w < (N-1)/N && continue

        λb = λ_b(g, w, N)
        lyap_vals_filled[j, i] = max(λa, λb, λ_d_realpart)
        #max(λa, λb, λ_d_realpart) > -1.0 && @info "largest eigenvalue at $(w),$(g): $(max(λa, λb,λ_d_realpart))"
    end
end

pl_lyaps_filled = plot_phase_diagram(lyap_vals_filled, ws_phasediag, gs_phasediag, N;
    xlims=(0.0, 1.0),
    ylims=(gs_phasediag[1], gs_phasediag[end]),
    critical_curve_colors=[:gray15, :gray15],
    critical_curve_lines=[:dash, :dot],
    critical_curve_lw = [3,3],
    HEATMAP_DEFAULTS...,
    ylabel=L"1/(T_x(N-1))", clims=(-1, 1),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    legend = false,
    colorbar = false,
    x_ticks=[0.0, 0.5, 1.0])


#---------------------read parameter values------------------

#horizontal direction on phasediagram 
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c_article/" #"c/"
params_orig_horizontal = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled_horizontal = readdlm(param_dir * "params_w_g.txt")
params_changed_horizontal = readdlm(param_dir * "params_c.txt")
cs = params_changed_horizontal
ws_horizontal = params_rescaled_horizontal[:,1]
gs_horizontal = params_rescaled_horizontal[:, 2]
g_horizontal = gs_horizontal[1]

#alpha (rp=rc) direction on phasediagram
param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_alpha_N.txt")
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = params_changed_curve[:,1]
rs_rounded = round.(rs;digits=3)
ws_curve_examples = params_rescaled_curve[:, 1]
gs_rcurve_examples = params_rescaled_curve[:, 2]


#------------------r curve--------------------------
r_curve(w, K, c) = K * c * (1 / w - 1) #g/(N-1) = 1/(ϵ*(N-1)), rp=rc
rc = 1.8;
K = 1.0;
α0 = 0.1;

_, w_B2, w_B4 = ws_horizontal
g_B2 = g_horizontal
c = g_B2 / K / (1 / w_B2 - 1)
r_curve(w_B2, K, c) #check if alpha curve passes through B3

w_B3,w_B2,w_B4 = ws_curve_examples
ws_B2_B4 = [w_B2-0.1:0.001:w_B4+0.1;]
gs_rcurve = r_curve.(ws_B2_B4, K, c) #alpha curve
g_B3, g_B2, g_B4 = gs_rcurve_examples

#plot rp=rc=r curve
plot!(ws_B2_B4, gs_rcurve, ls=:dot, la=1,
    lw=3, lc=:white, label="")

#------------------------plot param values for lifetimes---------------------

plot!([ws_horizontal[1]], [gs_horizontal[1]],
    ms=ms, st=:scatter, mc=:cyan,
    markerstrokewidth=msw,markershape=:diamond,
    label=L"c=%$(cs[1]), r=%$(rs[2])") 

annotate!(ws_horizontal[1] + 0.05, g_horizontal * 0.9, text(L"B_1", :cyan, :center, annot_size))

#------------------------plot param values for time series examples-------------------------
ws_example = [0.1, ws_horizontal[2], ws_horizontal[3], 0.95]
color_default = :gray15
example_colors = [:orangered :lime color_default :darkblue]
markershapes = [:circle :dtriangle :star :circle]

hline!([g_horizontal], ls=:dot, la=0.5, lw=2, lc=:gray10, label="") #horizontal line

plot!(ws_example', fill(g_horizontal, 4)',
    ms=ms, st=:scatter, mc=example_colors,
    markerstrokewidth=msw,
    markershape=markershapes, 
    label=["" L"c=%$(cs[2]),r=%$(rs_rounded[2])" L"c=%$(cs[3]),r=%$(rs_rounded[2])" ""],
    legend_columns=1,legendfontsize=8)

#add annotations for example time series
labels = [L"A", L"B_2", L"C", L"D"]
x_annot_offsets = [0.0, 0.05, 0.0, 0.0]
y_annot_offsets = [0.0, 1.4, 0.0, 0.0]
for (i, pl) in enumerate(ws_example)
    xpos = ws_example[i] + x_annot_offsets[i]
    ypos = g_horizontal * 0.9 + y_annot_offsets[i]
    annotate!(xpos, ypos, text(labels[i], example_colors[i], :center, annot_size))
end


plot!([w_B3 w_B4], [g_B3 g_B4],
    ms=ms, st=:scatter, mc=[:orange :purple],
    markershape=[:utriangle :rect],
    label=[L"c=%$(cs[2]),r=%$(rs_rounded[1])" L"c=%$(cs[2]), r=%$(rs_rounded[3])"])

annotate!(w_B3 + 0.08, g_B3 - 0.5, text(L"B_3", :orange, :center, annot_size))
annotate!(w_B4, g_B4 + 0.5, text(L"B_4", :purple, :left, annot_size))

#add colorbar
cbar_range = [-0.7:0.01:0.7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-0.5,0.0,0.5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-6Plots.mm,tickfontsize=tick_fontsize_large-5,
    top_margin=10Plots.mm)

#combine heatmap with the separate colorbar
l = @layout [a{0.95w} b]
pl_phase_diag = plot!(pl_lyaps_filled, pl_cbar, layout=l)

#add fig level annotation 
annotate!(pl_phase_diag,-0.15, 10.9, text(L"(c)", :gray10, :center, annot_size))

savefig(pl_phase_diag, figs_dir * "N10_phasediagram_extended5_80x80_filled_Ttr_10000_" * "labeled_horizontal_alpha_curve_article_params.png")
