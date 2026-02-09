using Plots, LaTeXStrings
using AttentionDynamics
using DelimitedFiles
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#-------------------------plot phase diagram---------------------

#--------------read lyaps--------------------
data_dir = "../../data/phasediag_data_N10_extended_d0_1.0e-9_ut1e-6/"
lyap_vals = reshape(readdlm(data_dir*"lyap_values.txt"),(60,60))
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[:,1]
gs = params[:,2]
N=10

pl_lyaps = plot_phase_diagram(lyap_vals,ws,gs,N;HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}",clims=(-0.5,0.5),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar_ticks = [-0.5,0.0,0.5])

#fill in with analytical eigenvalues
λa = -1
λ_b(g,w,N) = g*(1-1/(N-1)*1/(1/w-1)) 
lyap_vals_filled = deepcopy(lyap_vals)
for (j,g) in enumerate(gs) 
    for (i,w) in enumerate(ws)
        w < (N-1)/N && continue

        λb = λ_b(g, w, N)
        lyap_vals_filled[j, i] = max(λa, λb)
        max(λa, λb) > -1.0 && @info "largest eigenvalue at $(w),$(g): $(max(λa, λb))"
    end
end


pl_lyaps_filled = plot_phase_diagram(lyap_vals_filled, ws, gs, N; HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}", clims=(-0.5, 0.5),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large,
    colorbar_title="",
    colorbar_ticks=[-0.5, 0.0, 0.5])

ws_example = [0.5, 0.65, 0.75, 0.92]
g_example = 15.75
hline!([g_example / (N - 1)], ls=:dot, la=0.5, lw=2, lc=:gray10, label="")
plot!(ws_example, fill(g_example / (N - 1), 4),
    ms=5, st=:scatter, mc=:orange,markerstrokecolor=:orange,
    label="")

#add annotations 
labels = ["A", "B", "C", "D"]
x_annot_offsets = [0.01, 0.0, 0.0, -0.05]
for (i, pl) in enumerate(ws_example)
    xpos = ws_example[i] + x_annot_offsets[i]
    ypos = g_example/(N-1) * 1.05
    annotate!(pl_lyaps_filled, xpos, ypos, text(labels[i], :gray10, :center, 18))
end
pl_lyaps_filled

#add original param curves
αs = range(0.08,0.15;length=20)
rp = 2.0; rc = 2.0; K = 1.0; c = 3.0; α = 0.1
w_g_pairs_alpha = rescaled_params.(rp, rc, K, c, αs, N)
ws = [w_g_pairs_alpha[i][1] for i in 1:length(w_g_pairs_alpha)]
gs = [w_g_pairs_alpha[i][2] for i in 1:length(w_g_pairs_alpha)]

color = :lime
plot!(pl_lyaps_filled, w_c, g_c ./ (N - 1),
    mc=:lime, lc=:lime, lw=1, arrow=:arrow,ms=4,
    markershape=:circle,markerstrokecolor=:lime,
    label=L"\alpha \in [%$(αs[1]),%$(αs[end])]",
    legend=:topright)

savefig(pl_lyaps_filled,figs_dir*"N10_phasediagram_extended_filled_Ttr_1000"*"labeled.png")



#=
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