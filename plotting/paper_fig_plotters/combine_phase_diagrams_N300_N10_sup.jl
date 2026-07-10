include("plot_phase_diag_N300.jl") #N300
pl_phase_diag_N300 = pl_phase_diag #rename here 

#readd B1 point
plot!(pl_phase_diag_N300[1],[w_B1], [g_B2], ms=8, st=:scatter, mc=:cyan,label="",markershape=:diamond)
annotate!(pl_phase_diag_N300[1],w_B1 + 0.03, g_B2 + 0.6, text(L"B_1", :cyan, :center, annot_size_letters))


include("plot_phase_diag_N10.jl") #N10 
pl_phase_diag_N10 = pl_phase_diag #rename here 
plot!(pl_phase_diag_N10[1],[w_B3 w_B4], [g_B3 g_B4],
    ms=ms, st=:scatter, mc=[:orange :purple],
    markershape=[:utriangle :rect],label="")
#plot rp=rc=r curve
plot!(pl_phase_diag_N10[1],ws_B2_B4, gs_rcurve, ls=:dot, la=1,
    lw=3, lc=:white, label="")

plot!(pl_phase_diag_N10[1],[w_B3 w_B4], [g_B3 g_B4],
    ms=ms, st=:scatter, mc=[:orange :purple],
    markershape=[:utriangle :rect])

annotate!(pl_phase_diag_N10[1],w_B3 + 0.08, g_B3 - 0.5, text(L"B_3", :orange, :center, annot_size))
annotate!(pl_phase_diag_N10[1],w_B4, g_B4 + 0.5, text(L"B_4", :purple, :left, annot_size))

#add back B1
plot!([ws_horizontal[1]], [gs_horizontal[1]],
    ms=ms, st=:scatter, mc=:cyan,
    markerstrokewidth=msw,markershape=:diamond,
    label=L"c=%$(cs[1]), r=%$(rs[2])") 

annotate!(pl_phase_diag_N10[1],ws_horizontal[1] + 0.01, g_horizontal * 1.1, text(L"B_1", :cyan, :center, annot_size))
annotate!(pl_phase_diag_N10[1],w_B3 + 0.08, g_B3 - 0.5, text(L"B_3", :orange, :center, annot_size))
annotate!(pl_phase_diag_N10[1],w_B4, g_B4 + 0.5, text(L"B_4", :purple, :left, annot_size))


pl_phase_diag_N300_N10 = plot(pl_phase_diag_N10,pl_phase_diag_N300, 
    layout = (1,2), size = (1000,400),bottom_margin=6Plots.mm,
    top_margin=8Plots.mm,dpi=300)

#add fig level annotation 
annotate!(pl_phase_diag_N300_N10[1],-0.2, 10.8, text(L"(a)", :gray10, :center, annot_size))
annotate!(pl_phase_diag_N300_N10[3], -0.2, 10.8, text(L"(b)", :gray10, :center, annot_size))

savefig(pl_phase_diag_N300_N10,"../../paper_figs/phase_diagrams_N300_N10_sup.png")