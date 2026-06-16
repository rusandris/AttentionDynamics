include("plot_phase_diag_N300.jl")
include("plot_lyap_curves_largesys.jl")

#add inset to phase diag 

plot!(pl_phase_diag,ws_r, Λs_r, inset = (1, bbox(0.65,0.05,0.3,0.3)), subplot = 3, 
    lc=:gray10, ls=:dot, ylims=(2,9), ylabel=L"\Lambda",yticks=[3,6,9],xticks=[0.5,0.35,0.2],
        xlabel=L"w", guidefontsize=25,label="",lw=2,xlims=(0.15,0.5))

plot!(pl_phase_diag[3],[ws_r[28]], [Λs_r[28]], ms=6, st=:scatter, mc=:lime,markershape=:circle,label="")
#plot!(pl_phase_diag[3],[ws_r[end]], [Λs_r[end]], ms=6, st=:scatter, mc=:orange,markershape=:utriangle,label="")
plot!(pl_phase_diag[3],[ws_r[1]], [Λs_r[1]], ms=6, st=:scatter, mc=:purple,markershape=:rect,label="")


l = @layout [a{0.5w} b{0.5w}]
pl_fig2 = plot(pl_phase_diag,pl_lyaps, layout = l, size=(1000,400), 
    guidefontsize=25,tickfontsize=15,bottom_margin=5Plots.mm,
    top_margin=12Plots.mm,dpi=300)

plot!(pl_fig2[3],guidefontsize=20,tickfontsize=8) #finally adjust inset fontsizes
plot!(pl_fig2[4],right_margin=4Plots.mm) #adjust right margin
plot!(pl_fig2[2],right_margin=8Plots.mm) #adjust right margin


#add fig annotations
annotate!(pl_fig2[1], -0.14, 11.5, text("(a)", :gray10, :right, annot_size))
annotate!(pl_fig2[4], -0.1, 9.2, text("(b)", :gray10, :right, annot_size))


savefig(pl_fig2,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_lyap_curves_inset.png")