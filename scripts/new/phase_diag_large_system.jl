using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("./bif_floq_od_method_parameters.jl")



#read w and g params

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-2,adaptive=true,maxiters=1e15)
N = 300
ws = range(0.1,0.65;length=50)
g_rel = 5
gs = range(4.0,6.0;length=50) .* (N-1)


#try original equations
u0 = 0.1 .* rand(2N)
rc=9.0;rp=9.0;K=1.0;c=2.4;α=0.005
p = [rc,rp,K,c,α]
ds = create_original_ds(u0,p)

ws_crit = (N-1) ./ gs
ws = [0.01:0.01:0.99;]
pl = plot( ws,1 ./ ws , lc=:gray10,lw=1, ylabel=L"g/(N-1)",xlabel=L"w", 
    guidefontsize=15,label="")

plot!()



wc2 = (N-1)/N

#generate standard init cond
T_fp = 100
Ttr_fp = 100
Δt_fp = 0.1
fp_settings = [T_fp,Ttr_fp,Δt_fp]
Ttr = 1000
Ttrs = [Ttr]
T = 5000
lyap_vals = []


@time for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds,gs,ws; T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
        d0=1e-9, ut=1e-6, lt=1e-12,save=false,realtime_plotting=false,parallelize=false)
    push!(lyap_vals,Λ)
end

writedlm(data_dir*"phasediag_params.txt",hcat(ws,gs))
writedlm(data_dir*"lyap_values.txt",lyap_vals)
