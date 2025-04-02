using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots
cd(@__DIR__)

#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients

#read w and g params
ws_saved = readdlm("../../plotting/data/N=10.0_wvalues.txt")
gs_saved = readdlm("../../plotting/data/N=10.0_gvalues.txt")
gs_N = gs_saved / (10-1)

diffeq = (alg = Vern9(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
N = 7
ds = create_ds(rand(2N),[N+1,0.5,N];diffeq=diffeq)
gs = gs_N * (N-1)
Ttrs = [3000,50000,100000]
T = 5000
g_val_test = gs[24] #calc for only this test values (one row)
set_parameter!(ds,1,g_val_test)
ws = ws_saved[ws_saved .> (N-1)/g_val_test]
lyap_vals = zeros(length(ws),length(Ttrs))

for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    Λ = lyapunov_diagram(ds,[g_val_test],ws; T=T, Ttr=Ttr, Δt=0.1, 
        d0=1e-9, ut=1e-6, lt=1e-12,save=false,realtime_plotting=false,parallelize=false)
    @show size(Λ)
    lyap_vals[:,i] .= Λ'  
end

#plot_phase_diagram(Λ_N6,ws_saved[:],gs_N6[:],5;HEATMAP_DEFAULTS...)
plot(ws,lyap_vals,ylims=(-0.2,Inf),labels = ["Ttrs=$(Ttrs[1])" "Ttrs=$(Ttrs[2])" "Ttrs=$(Ttrs[3])"])