using Plots, LaTeXStrings
using AttentionDynamics
using DelimitedFiles
cd(@__DIR__)
include("defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#plot settings 
annot_size = 20
ms = 6

#-------------------------plot phase diagram---------------------

#--------------read lyaps--------------------
figs_dir = "../../paper_figs/"
#data_dir = "../../data/phasediag_data_N10_extended_d0_1.0e-9_ut1e-6/"
#data_dir = "../../data/phasediag_data_N10_extended3_80x80/"
data_dir = "../../../../../Dropbox/att_dyn/data/data_phasediag/phasediag_data_N10_extended5_80x80/"
#lyap_vals = reshape(readdlm(data_dir*"lyap_values_Ttr10000.txt"),(80,80))
lyap_vals = reshape(readdlm(data_dir * "lyap_values_Ttr10000.txt"), (80, 80))
params = readdlm(data_dir*"phasediag_params.txt")
ws_phasediag = params[:,1]
gs_phasediag = params[:, 2]
N=10

pl_lyaps = plot_phase_diagram(lyap_vals, ws_phasediag, gs_phasediag, N;
    xlims=(0.0, 1.0),
    ylims=(gs_phasediag[1], gs_phasediag[end]) ./ (N - 1),
    critical_curve_colors=[:gray10,:gray10],
    critical_curve_lines=[:dash, :dot],
    HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}",clims=(-1,1),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    plot_critical_curve=true,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar = false,
    x_ticks = [0.0, 0.5, 1.0])


#fill in with analytical eigenvalues
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
    ylims=(gs_phasediag[1], gs_phasediag[end]) ./ (N - 1),
    critical_curve_colors=[:gray15, :gray15],
    critical_curve_lines=[:dash, :dot],
    critical_curve_lw = [3,3],
    HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}", clims=(-1, 1),
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
params_orig_curve = readdlm(param_dir * "params_K_c_N.txt")
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = params_changed_curve[:,1]
rs_rounded = round.(rs;digits=3)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2]

#------------------------plot param values for lifetimes---------------------

plot!([ws_horizontal[1]], [gs_horizontal[1]],
    ms=ms, st=:scatter, mc=:cyan,
    markerstrokewidth=0.0,markershape=:dtriangle,
    label=L"c=%$(cs[1]), r=%$(rs[2])") 

annotate!(ws_horizontal[1] + 0.05, g_horizontal * 0.9, text(L"B_1", :cyan, :center, annot_size))

#------------------------plot param values for time series examples-------------------------
ws_example = [0.1, ws_horizontal[2], ws_horizontal[3], 0.95]
color_default = :gray15
example_colors = [:orange :lime color_default :gray80]
markershapes = [:rect :circle :utriangle :diamond]

hline!([g_horizontal], ls=:dot, la=0.5, lw=2, lc=:gray10, label="") #horizontal line

plot!(ws_example', fill(g_horizontal, 4)',
    ms=ms, st=:scatter, mc=example_colors,
    markerstrokewidth=0.0, markerstrokecolor=example_colors,
    markershape=markershapes, 
    label=["" L"c=%$(cs[2]),r=%$(rs_rounded[2])" L"c=%$(cs[3]),r=%$(rs_rounded[2])" ""],
    legend_columns=1,legendfontsize=8)

#add annotations for example time series
labels = [L"A", L"B_2", L"C", L"D"]
x_annot_offsets = [0.0, 0.08, 0.0, 0.0]
y_annot_offsets = [0.0, 0.0, 0.0, 0.0]
for (i, pl) in enumerate(ws_example)
    xpos = ws_example[i] + x_annot_offsets[i]
    ypos = g_horizontal * 0.9 + y_annot_offsets[i]
    annotate!(xpos, ypos, text(labels[i], example_colors[i], :center, annot_size))
end

#------------------alpha curve--------------------------
alpha_curve(w, K, c) = K * c * (1 / w - 1) #g/(N-1) = 1/(ϵ*(N-1)), rp=rc
rc = 1.8;
K = 1.0;
α0 = 0.1;

_, w_B2, w_B3 = ws_horizontal
g_B2 = g_horizontal
c = g_B2 / K / (1 / w_B2 - 1)
alpha_curve(w_B2, K, c) #check if alpha curve passes through B3

w_B4 = ws_curve[1]
w_B5 = ws_curve[3]
ws_B4_B5 = [w_B4:0.001:w_B5;]
gs_alpha = alpha_curve.(ws_B4_B5, K, c) #alpha curve
g_B4 = gs_alpha[1]
g_B5 = gs_alpha[end]

#plot rp=rc=r curve
plot!(ws_B4_B5, gs_alpha, ls=:dot, la=1,
    lw=2, lc=:lightblue, label="")

plot!([w_B4 w_B5], [g_B4 g_B5],
    ms=ms + 4, st=:scatter, mc=[:blue :darkgreen],
    markerstrokewidth=0.0, markershape=[:ltriangle :rtriangle],
    label=[L"c=%$(cs[2]),r=%$(rs_rounded[1])" L"c=%$(cs[2]), r=%$(rs_rounded[3])"])

annotate!(w_B4 + 0.08, g_B4 - 0.5, text(L"B_3", :blue, :center, annot_size))
annotate!(w_B5, g_B5 + 0.5, text(L"B_4", :darkgreen, :left, annot_size))

#add colorbar
cbar_range = [-0.7:0.01:0.7;]
pl_cbar = heatmap([1], cbar_range, [cbar_range;;],yticks = [-0.5,0.0,0.5], 
    c=HEATMAP_DEFAULTS.c, fg_color=:gray10, colorbar=false, xaxis=false, 
    ymirror=true, left_margin=-6Plots.mm,tickfontsize=tick_fontsize_large-5)

#conbine heatmap with the separate colorbar
l = @layout [a{0.95w} b]
pl = plot!(pl_lyaps_filled, pl_cbar, layout=l, tick_direction=:out)

savefig(pl, figs_dir * "N10_phasediagram_extended5_80x80_filled_Ttr_10000_" * "labeled_horizontal_alpha_curve_article_params.png")


#=
#example parameters 

#ws_example = [0.3, 0.48, 0.7, 0.95]
#ws_example = [0.3, 0.63, 0.8, 0.95]
w_B2 = 0.48
ws_example = [0.3, w_B2, 0.7, 0.95]
g_example = 22.5 #9*1.75

=#
#=

=#


#=
plot!([w_B1,w_B3],[g_example, g_example] ./ (N-1),
    ms=ms, st=:scatter, mc=[:lime,:cyan],
    markerstrokewidth=0.0,
    label="")

annotate!(pl_lyaps_filled, w_B1, g_example/(N-1)*1.05, text(L"B_1", :lime, :center, 14))
annotate!(pl_lyaps_filled, w_B3-0.02, g_example/(N-1)*1.05, text(L"B_3", :cyan, :center, 14))
=#




#=
#----------------------------add original param curves---------------------------

#---------------rp------------------

rc = 1.8;
K = 1.0;
c0 = 2.5;
α0 = 0.1;


#B is one corner of the triangle
w_B2 = ws_example[2]
g_B2 = g_example
#rp_B1 = g_example*α0

#intersection of vertical and alpha curve 
g_alpha(w,rp,rc,K,c) = rp*K*c/rc * (1/w-1) #alpha curve

#check if alpha curve passes through B3
rp_B2 = rc/(K*c0)*g_B2/(1/w_B2-1)

g_B2 = g_alpha(w_B2, rp_B2, rc, K, c0)
w_B1 = 0.43
w_B3 = 0.53
w_B1_B3 = [w_B1:0.001:w_B3;]
gs_alpha = g_alpha.(w_B1_B3, rp_B2, rc, K, c0) #alpha curve
g_B1 = gs_alpha[1]
g_B3 = gs_alpha[end]
#plot triangle 
plot!([w_B1], [g_B3 / (N - 1)], ms=ms, mc=:lime, markerstrokewidth=0.0, st=:scatter, label="")
plot!([w_B1], [g_B2 / (N - 1)], ms=ms, mc=:cyan, markerstrokewidth=0.0, st=:scatter, label="")
plot!(w_B1_B3, gs_alpha ./ (N - 1), st=:scatter, ms=0.1, la=0.5, lw=2, lc=:gray10, label="")
plot!([w_B1, w_B1], [g_B3, g_B2] ./ (N - 1), ls=:dot, la=0.5, lw=2, lc=:gray10, label="")

annotate!(pl_lyaps_filled, w_B1 - 0.01, g_B3 / (N - 1) + 0.1, text(L"B_1", :lime, :right, 14))
annotate!(pl_lyaps_filled, w_B1, g_B2 / (N - 1) + 0.1, text(L"B_2", :cyan, :right, 14))


#save three B parameter pairs into file 
params_savedir = "../../data/lifetime_distributions_rescaled_alpha_curve/"
mkpath(params_savedir)
writedlm("../../data/lifetime_distributions_rescaled_alpha_curve/B_param_pairs.txt", hcat([w_B1, w_B2, w_B3], [g_B1, g_B2, g_B3]))

#plot!([w_B3], [g_B3 / (N - 1)], ms=6, mc=:cyan, st=:scatter)

savefig(pl_lyaps_filled, figs_dir * "N10_phasediagram_extended3_80x80_filled_Ttr_10000_" * "labeled_triangle.png")
=#

#=
rs = range(rp_B1, rp_B2; length=3)
w_g_pairs_rp = rescaled_params.(rs, rc, K, c0, α0, N)
ws_rp = [w_g_pairs_rp[i][1] for i in 1:length(w_g_pairs_rp)]
gs_rp = [w_g_pairs_rp[i][2] for i in 1:length(w_g_pairs_rp)]


#---------------alpha------------------
#only coincides with r curve if rp=rc
αs = range(α_B1, α_B3; length=3)
w_g_pairs_alpha = rescaled_params.(rp_B2, rc, K, c0, αs, N)
ws_alpha = [w_g_pairs_alpha[i][1] for i in 1:length(w_g_pairs_alpha)]
gs_alpha = [w_g_pairs_alpha[i][2] for i in 1:length(w_g_pairs_alpha)]

#---------------c------------------
c_B3 = c0
c_B1 = c0

cs = range(3.0, 3.0/αs[end]*αs[1]; length=5)
rp = 2.0;rc = 2.0; K = 1.0; α = αs[end];
w_g_pairs_cs = rescaled_params.(rp, rc, K, cs, α, N)
ws_c = [w_g_pairs_cs[i][1] for i in 1:length(w_g_pairs_cs)]
gs_c = [w_g_pairs_cs[i][2] for i in 1:length(w_g_pairs_cs)]

#plot triangle 
plot!(pl_lyaps_filled, ws_rp, gs_rp ./ (N - 1),
    mc=[:darkgreen, :darkgreen, :lime],
    markershape=:circle,lc=:darkgreen,
    markerstrokewidth=0.0,
    ms=[6,3,6],label="")

plot!(pl_lyaps_filled, ws_alpha, gs_alpha ./ (N - 1),
    mc=[:lime, :lime, :lime, :lime, :cyan],
    markershape=:circle, lc=:lime,
    markerstrokewidth=0.0,
    ms=[6,3,6], label="")

plot!(pl_lyaps_filled, ws_c, gs_c ./ (N - 1),
    mc=[:cyan, :cyan, :cyan, :cyan, :darkgreen],
    markershape=:circle, lc=:cyan,
    markerstrokewidth=0.0,
    ms=[6,3,6], label="")



#---------------r------------------

rs = range(1.0, 3.0; length=20)
K = 1.0; c = 3.0;α = 0.1;
w_g_pairs_rs = rescaled_params.(rs, rs, K, c, α, N)
ws = [w_g_pairs_rs[i][1] for i in 1:length(w_g_pairs_rs)]
gs = [w_g_pairs_rs[i][2] for i in 1:length(w_g_pairs_rs)]

color = :orange
plot!(pl_lyaps_filled, ws, gs ./ (N - 1),
    mc=color, lc=color, lw=1, arrow=:arrow, ms=4,
    markershape=:circle, markerstrokecolor=color,
    label=L"r_p = r_c \in [%$(rs[1]),%$(rs[end])]",
    legend=:topright)

#---------------rc------------------

rs = range(1.0, 3.0; length=20)
rp = 2.0; K = 1.0;c = 3.0;α = 0.1;
w_g_pairs_rc = rescaled_params.(rp, rs, K, c, α, N)
ws = [w_g_pairs_rc[i][1] for i in 1:length(w_g_pairs_rc)]
gs = [w_g_pairs_rc[i][2] for i in 1:length(w_g_pairs_rc)]

color = :blue3
plot!(pl_lyaps_filled, ws, gs ./ (N - 1),
    mc=color, lc=color, lw=1, arrow=:arrow, ms=4,
    markershape=:circle, markerstrokecolor=color,
    label=L"r_p = %$(rp),r_c \in [%$(rs[1]),%$(rs[end])]",
    legend=:topright)

plot!(pl_lyaps_filled, legendfontsize=10)


savefig(pl_lyaps_filled,figs_dir*"N10_phasediagram_extended_filled_Ttr_1000"*"labeled.png")

#-------------------plot original parameter curves on phase diagram-------------

N=300
rc = 9.0;
rp = 9.0;
K = 1.0;
c = 1.4;
α = 0.005;
α2 = 0.005

#plot critical curve
ws = [0.01:0.01:0.99;]
pl_crit_curve = plot(ws, 1 ./ ws, lc=:gray10, lw=1, 
    ylabel=L"g/(N-1)", xlabel=L"w",
    guidefontsize=15, label=L"g(w)^H = 1/w ",
    title = L"N=%$(N), K=%$(K), c=%$(c)",
    legendfontsize=12)

r_test = 6.0
#rs = [6.0:10.0:40.0;]
#rs = [6.0:1.0:12.0;]
rs = [6.0:3.0:16.0;]


#when rp=rc
w_g_rs = rescaled_params.(rs, rs, K, c, α, N)
w_rs = [w_g_rs[i][1] for i in 1:length(w_g_rs)]
g_rs = [w_g_rs[i][2] for i in 1:length(w_g_rs)]

#plot rp = rc = r params 
plot!(pl_crit_curve, w_rs, g_rs ./ (N - 1),
    mc=:green, lc=:green, lw=1, arrow=:arrow,
    label=L"r_p=r_c \in [%$(rs[1]),%$(rs[end])]", ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)
#when rp \neq rc, only rp changes
w_g_rp = rescaled_params.(rs, r_test, K, c, α, N)
w_rp = [w_g_rp[i][1] for i in 1:length(w_g_rp)]
g_rp = [w_g_rp[i][2] for i in 1:length(w_g_rp)]

plot!(pl_crit_curve, w_rp, g_rp ./ (N - 1),
    mc=:purple, lc=:purple, lw=1, arrow=:arrow,
    label=L"r_p \in [%$(rs[1]),%$(rs[end])]",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

#when rp \neq rc, only rc changes
rs = [4.0:0.5:6.0;]
w_g_rc = rescaled_params.(6.0, rs, K, c, α, N)
w_rc = [w_g_rc[i][1] for i in 1:length(w_g_rc)]
g_rc = [w_g_rc[i][2] for i in 1:length(w_g_rc)]

#plot rp = rc = r params 
plot!(pl_crit_curve, w_rc, g_rc ./ (N - 1),
    mc=:blue, lc=:blue ,lw=1, arrow=:arrow,
    label=label=L"r_c \in [%$(rs[1]),%$(rs[end])]",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

ylims!(pl_crit_curve,(0,10))
xlims!(pl_crit_curve, (0.1, 0.4))
savefig(pl_crit_curve,"../../paper_figs/")

#change cs
r_test = 16.0
cs = [1.5:0.3:2.5;]

w_g_c = rescaled_params.(r_test, r_test, K, cs, α, N)
w_c = [w_g_c[i][1] for i in 1:length(w_g_c)]
g_c = [w_g_c[i][2] for i in 1:length(w_g_c)]

plot!(pl_crit_curve, w_c, g_c ./ (N - 1),
    mc=:blue, lc=:blue ,lw=1, arrow=:arrow,
    label=label=L"c \in [%$(cs[1]),%$(cs[end])]",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

=#

#=
#add more points to B phase

#w direction
w_B1 = 0.43
w_B3 = 0.53
g_B4 = g_example

plot!([w_B1, w_B3], [g_B4, g_B4] ./ (N - 1),
    ms=ms, st=:scatter, mc=[:lime, :cyan],
    markerstrokewidth=0.0,
    label="")

annotate!(pl_lyaps_filled, w_B4, g_B4 / (N - 1) + 0.1, text(L"B_1", :lime, 14))
annotate!(pl_lyaps_filled, w_B5 - 0.01, g_B4 / (N - 1) + 0.1, text(L"B_3", :cyan, 14))

#g direction
g_B4 = g_example - 1.5
g_B5 = g_example + 1.5

plot!([w_B2, w_B2], [g_B4, g_B5] ./ (N - 1),
    ms=ms, st=:scatter, mc=[:purple, :blue],
    markerstrokewidth=0.0,
    label="")

plot!([w_B2, w_B2], [g_B4, g_B5] ./ (N - 1), ls=:dot, lc=:gray10, label="")
annotate!(pl_lyaps_filled, w_B2 + 0.03, g_B4 / (N - 1) - 0.1, text(L"B_4", :purple, :center, 14))
annotate!(pl_lyaps_filled, w_B2 - 0.02, g_B5 / (N - 1) + 0.1, text(L"B_5", :blue, :center, 14))


#save param pairs into file 
savedir_param = "../../data/trajs_precomp_rescaled_horizontal_vertical/"
writedlm(savedir_param * "param_pairs.txt", [g_example w_B1; g_example w_B2; g_example w_B3; g_B4 w_B2; g_B5 w_B2])
=#