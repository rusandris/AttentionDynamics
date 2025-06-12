using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")


#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients

#read w and g params
#ws_saved = readdlm("../../plotting/data/N=10.0_wvalues.txt")
#gs_saved = readdlm("../../plotting/data/N=10.0_gvalues.txt")
#gs = gs_saved / (10-1)

data_dir = "../../phasediag_data_300/"
mkpath(data_dir)
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
N = 300
ws = range(0.4,0.75;length=20)
gs = range(1800,2000;length=20)

wc1 = (N-1)/gs[end]
wc2 = (N-1)/N

#generate standard init cond
#u0 = rand(2N)
u0 = generate_initial_condition(N)
ds = create_ds(u0,[gs[1],0.285,N];diffeq=diffeq)
T_fp = 100
Ttr_fp = 100
Δt_fp = 0.1
fp_settings = [T_fp,Ttr_fp,Δt_fp]
Ttrs = [1000]
T = 5000
lyap_vals = []

@time for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds,gs,ws; T=T, Ttr=Ttr, Δt=1.0, fp_settings=fp_settings, 
        d0=1e-9, ut=1e-6, lt=1e-12,save=false,realtime_plotting=false,parallelize=false)
    push!(lyap_vals,Λ)
end

for (i,Λ) in enumerate(lyap_vals)
    writedlm("../../N5_phasediagram_test_Ttr_max_$(Ttrs[i])"*".txt",Λ)

end
writedlm("../../N5_phasediag_params.txt",hcat(ws,gs))
pl1 = plot_phase_diagram(lyap_vals[1],ws[:],gs[:],5;HEATMAP_DEFAULTS...)
pl2 = plot_phase_diagram(lyap_vals[2],ws_saved[:],gs_N[:],5;HEATMAP_DEFAULTS...)
pl3 = plot_phase_diagram(lyap_vals[3],ws_saved[:],gs_N[:],5;HEATMAP_DEFAULTS...)
savefig(pl1,"../../N5_phasediagram_Ttr_max_$(Ttrs[1])"*".png")
savefig(pl2,"../../N5_phasediagram_Ttr_max_$(Ttrs[2])"*".png")
savefig(pl3,"../../N5_phasediagram_Ttr_max_$(Ttrs[3])"*".png")


plot(ws,lyap_vals,ylims=(-0.2,Inf),labels = ["Ttrs=$(Ttrs[1])" "Ttrs=$(Ttrs[2])" "Ttrs=$(Ttrs[3])"])