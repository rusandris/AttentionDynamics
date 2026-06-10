include("plot_phase_diag_N300.jl")
include("plot_lyap_curves_largesys.jl")

#add inset to phase diag 

plot!(pl_pd300,ws_r, Λs_r, inset = (1, bbox(0.65,0.05,0.3,0.3)), subplot = 2, 
    lc=:gray10, ylims=(2,9), xflip=true, ylabel=L"\Lambda",yticks=[3,6,9],xticks=[0.5,0.35,0.2],
        xlabel=L"w", guidefontsize=25,label="",lw=2,xlims=(0.15,0.5))

plot!(pl_pd300[2],[ws_r[28]], [Λs_r[28]], ms=6, st=:scatter, mc=:lime,markershape=:circle,label="")
plot!(pl_pd300[2],[ws_r[end]], [Λs_r[end]], ms=6, st=:scatter, mc=:orange,markershape=:utriangle,label="")
plot!(pl_pd300[2],[ws_r[1]], [Λs_r[1]], ms=6, st=:scatter, mc=:purple,markershape=:rect,label="")


#layout = (1,2)
pl_fig2 = plot(pl_pd300,pl_lyaps, layput = (1,2), size=(1200,500), 
    guidefontsize=25,tickfontsize=15,bottom_margin=5Plots.mm,left_margin=8Plots.mm,right_margin=5Plots.mm,dpi=300)


savefig(pl_fig2,"../../paper_figs/phase_diag_N300_Ttrmax_500_times_scaled_unfinished_lyap_curves_inset.png")