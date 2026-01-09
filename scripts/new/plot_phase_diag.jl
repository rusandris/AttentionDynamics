using Plots, LaTeXStrings
using AttentionDynamics

#-------------------------plot phase diagram---------------------
#check lifetime distributions while keeping rp/alpha constant
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

#TODO: 
#map different parts of the phase diagram
#plot lyapunov exp, avg \tau and avg. inequality