using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("../../scripts/new/bif_floq_od_method_parameters.jl")


#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients

#read w and g params
#ws_saved = readdlm("../../plotting/data/N=10.0_wvalues.txt")
#gs_saved = readdlm("../../plotting/data/N=10.0_gvalues.txt")
#gs = gs_saved / (10-1)


#data_dir = "./phasediag_data_N10/"
data_dir = "./phasediag_data_N10_zoomed/"
mkpath(data_dir)
figs_dir = "./poster_figs/"
mkpath(figs_dir)

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-2,adaptive=true,maxiters=1e15)
N = 10
ws = range(0.5,0.85;length=50)
g_rel = 5
gs = range(1.2,3.0;length=50) .* (N-1)
#=
#try original equations
u0 = generate_initial_condition(N)
rc=9.0;rp=9.0;K=1.0;c=2.4;α=0.005
p = [rc,rp,K,c,α]
ds = create_original_ds(u0,p)
=#
wc1 = (N-1)/gs[1]
wc2 = (N-1)/N

#generate standard init cond
#u0 = rand(2N)
u0 = generate_initial_condition(N)
ds = create_ds(u0,[5,0.5,N];diffeq)
T_fp = 100
Ttr_fp = 100
Δt_fp = 0.1
fp_settings = [T_fp,Ttr_fp,Δt_fp]
Ttrs = [Ttr]
T = 5000
lyap_vals = []

#
@time for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds,gs,ws; T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
        d0=1e-9, ut=1e-6, lt=1e-12,save=false,realtime_plotting=false,parallelize=false)
    push!(lyap_vals,Λ)
end

writedlm(data_dir*"phasediag_params.txt",hcat(ws,gs))
writedlm(data_dir*"lyap_values.txt",lyap_vals)


#=
#--------------read lyaps--------------------
lyap_vals = reshape(readdlm(data_dir*"lyap_values.txt"),(40,40))
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[:,1]
gs = params[:,2]
N=10
=#

pl1 = plot_phase_diagram(lyap_vals[1],ws,gs,N;HEATMAP_DEFAULTS...,
    ylabel=L"1/(\epsilon(N-1))",right_margin=3Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xguidefontsize=guide_fontsize_large)

#read in NS values
readdir("../../")
ns_bifpoints_data = readdlm("../../data/NS_bifpoints_data_gN10/contres_po_bifp_ns.txt")
plot!(pl1,ns_bifpoints_data[:,2],gs[1:6] ./ 9,lc=:gray10,mc=:gray10,ms=2)
#pl2 = plot_phase_diagram(lyap_vals[2],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
#pl3 = plot_phase_diagram(lyap_vals[3],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
savefig(pl1,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])"*".png")
savefig(pl2,"../../N5_phasediagram_Ttr_max_$(Ttrs[2])"*".png")
savefig(pl3,"../../N5_phasediagram_Ttr_max_$(Ttrs[3])"*".png")


plot(ws,lyap_vals,ylims=(-0.2,Inf),labels = ["Ttrs=$(Ttrs[1])" "Ttrs=$(Ttrs[2])" "Ttrs=$(Ttrs[3])"])
