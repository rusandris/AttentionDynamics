include("plot_phase_diag_N.jl")
include("plot_lyapunovs_effectsize.jl")

l = @layout [a{0.65w} b{0.35w}]
pl_pd_chaos_ends = plot(pl_phasediag_N,pl_lyap_curves_chaos_ends,layout=l,size=(800,400),
    guidefontsize = 22,tickfontsize = 12,bottom_margin=5Plots.mm)

#add fig level annotations
annotate!(pl_pd_chaos_ends[1], -0.2, 11.2, text(L"(a)", :gray10, :center, annot_size))
annotate!(pl_pd_chaos_ends[2], 0.02, 4.3, text(L"(b)", :gray10, :center, annot_size))
annotate!(pl_pd_chaos_ends[2], 0.02, -2.0, text(L"(c)", :gray10, :center, annot_size))



savefig(pl_pd_chaos_ends,"../../paper_figs/supplimentary_system_size_effect_with_chaos_endpoints_T1e3_Ttr5x1e4.png")
