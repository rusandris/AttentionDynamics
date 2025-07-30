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