using Revise
using DynamicalSystems
using DelimitedFiles
using AttentionDynamics
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")

#data_dir = "../../data/phasediag_data_N10/"
data_dir = "../../data/phasediag_data_N10_zoomed/"
#mkpath(data_dir)
figs_dir = "../../paper_figs/"
#mkpath(figs_dir)


#-----------------------------read lyaps-------------------------------
lyap_data = collect(readdlm(data_dir*"lyap_values.txt"))
n = Int(sqrt.(length(lyap_data)))
lyap_vals = reshape(readdlm(data_dir*"lyap_values.txt"),(n,n))
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[:,1]
gs = params[:,2]
N=10

pl1 = plot_phase_diagram(lyap_vals,ws,gs,N;HEATMAP_DEFAULTS...,
    ylabel=L"1/(\epsilon(N-1))",right_margin=3Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xguidefontsize=guide_fontsize_large,
    legend=:topright)

#plot original parameter spaces

rp = 1.0;rc = 1.0;K = 1.0;

#plot (g(α),w(α)) curve 

αs = [6*1e-2:1e-3:1e-1;]
ws_crit1_alpha = (N-1) ./ ( rp ./ αs )

w_g_alpha = rescaled_params.(rp,rc,c,K,αs,N)
w_g_alpha = [[w_g_alpha[i][1],w_g_alpha[i][2]] for i in 1:length(αs)]
w_alpha = [w_g_alpha[i][1] for i in 1:length(w_g_alpha)]
g_alpha = [w_g_alpha[i][2] for i in 1:length(w_g_alpha)]

#plot difference w - w_crit if alpha is changed
#plot(αs,(w_alpha - ws_crit1_alpha),lc=:purple,lw=3,ylabel=L"w_c(\alpha)  - w(\alpha)",
#    xlabel=L"\alpha",guidefontsize=15)

#plot distance to criticality (g_crit,w_crit) if alpha is changed
gs_crit1_alpha = (N-1) ./ ws_crit1_alpha
w_g_alpha_crit = [ [ws_crit1_alpha[i],gs_crit1_alpha[i]] for i in 1:length(ws_crit1_alpha)]
dist_crit = norm.( w_g_alpha .- w_g_alpha_crit  )
plot(αs,dist_crit,lc=:purple,lw=3,ylabel="distance to criticality",
    xlabel=L"\alpha",guidefontsize=15)

#plot (g(c),w(c)) curve
cs = [2.5:0.05:4.0;]
w_g_cs = rescaled_params.(rp,rc,cs,K,αs[1],N)
w_cs = [w_g_cs[i][1] for i in 1:length(w_g_cs)]
g_cs = [w_g_cs[i][2] for i in 1:length(w_g_cs)]

plot!(pl1,w_cs,g_cs ./ (N-1),lc=:green,lw=3,arrow=:arrow)







#=
#cs = [2.5:1.0:7.0;]
cs = [3.0,4.0]
linecolors = [:cyan,:purple]
#linealphas = [0.2:0.2:1.0;]
markershapes = [:circle,:rect,:cross,:x,:utriangle]
for (i,c) in enumerate(cs)
    curve = w_g_dependence_curve_adim.(ws;c=c)
    plot!(pl1,ws,curve;lw=1,
        lc=linecolors[i],label=L"g(w)_{c=%$c}",
        mc=:green)
end
pl1
#savefig(pl1,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])"*".png")


#-------------plot parameters for original parameter interpretation------------

#original params
N = 10
rp = 1.0;rc = 1.0;K = 1.0;α = 5*1e-2
c1 = 3.0;c2 = 4.0

w1,g1 = rescaled_params(rp,rc,K,c1,α,N)
w2,g2 = rescaled_params(rp,rc,K,c2,α,N)

plot!([w1],[g1/(N-1)],st=:scatter,lc=:cyan,ms=6,label="")
plot!([w2],[g2/(N-1)],st=:scatter,lc=:purple,ms=6,label="")

figs_dir = "../../paper_figs/"
savefig(pl1,figs_dir*"phasediag_data_N10_zoomed_orig_param_markers.png")
=#