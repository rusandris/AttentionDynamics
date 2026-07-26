using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
#using Plots,LaTeXStrings
cd(@__DIR__)

#--------------------system setup-------------------------
N = 10
g = 45.0
w = 0.4
p = [g,w,N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
u0 = rand(2N)
ds = create_ds(u0,p;diffeq=diffeq)

#--------------------ensemble-------------------------
n_ens = 1000
init_states = [rand(2N) for i in 1:n_ens]
Tmax = 1000000 #max lifetime
#--------------------list of parameters----------
w_NS = 0.329206866618409 #Neimark-Sacker for N=10, g=45 (1/eps(N-1) = 5.0)
#ws = [0.4:-0.001:w_NS+0.001;]
#ws = range(0.4,0.33,length=60)
ws_max = 0.4
Δws_min = 0.001 
Δws_max = ws_max - w_NS
Δws = logrange(Δws_min,Δws_max;length=60)
ws = reverse(w_NS .+ Δws)

#ws = reverse(logrange(w_NS+Δws_min,0.4;length=60))
#----------------------call lifetime calc on ensemble------------------
savedir = "../../../../../Insync/rusandris@gmail.com/Dropbox/att_dyn/data/data_escape_times/transient_lifetimes_from_ensembles_logdiffparams/"
mkpath(savedir)

writedlm(savedir*"w_params.txt", ws)
writedlm(savedir*"w_params_diffs.txt", Δws)
writedlm(savedir*"ens_size_Tmax.txt",(n_ens,Tmax))

@time for (i,w) in enumerate(ws) 
    @show (i,w)
    set_parameter!(ds,2,w)
    @time lts = transient_lifetimes(ds,init_states;Tmax=Tmax)
    string_zeros = i ÷ 10 == 0 ? "0" : ""
    writedlm(savedir*"lifetimes_$(string_zeros)$(i)_w$(w).txt",lts)
end

