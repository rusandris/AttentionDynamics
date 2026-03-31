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
n_ens = 10
init_states = [rand(2N) for i in 1:n_ens]

#--------------------list of parameters----------
w_NS = 0.329206866618409 #Neimark-Sacker for N=10, g=45 (1/eps(N-1) = 5.0)
ws = [0.4:-0.001:0.35;]

#----------------------call lifetime calc on ensemble------------------
savedir = "../../../data/transient_lifetimes_from_ensembles_test/"
mkpath(savedir)

for (i,w) in enumerate(ws) 
    @show w
    set_parameter!(ds,2,w)
    lts = transient_lifetimes(ds,init_states;Tmax=10000)
    string_zeros = i ÷ 10 == 0 ? "0" : ""
    writedlm(savedir*"lifetimes_$(string_zeros)$(i)_w$(w).txt",lts)
end

