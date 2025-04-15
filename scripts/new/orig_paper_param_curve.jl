using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/plot_functions/lyapunov.jl")

function dist_to_fixedpoint(u,fp)
    return norm(u - fp)
end

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
ws = [0.3:0.01:0.9;]
gs = w_g_dependence_curve_adim.(ws)

Ns = [3,4,5,7,10]
#N = 5

T = 2000
Ttr = 1000
lyaps_Ns = []

for N in Ns
    @show N
    ds = create_ds(rand(2N),[gs[1],ws[1],N];diffeq=diffeq)
    lyaps = zeros(length(ws))
    for i in 1:length(ws)
        @show ws[i]
        set_parameter!(ds,1,gs[i]*(N-1))
        set_parameter!(ds,2,ws[i])
        lyaps[i] = lyapunov(ds,T;Ttr=Ttr,Δt=0.5)
    end
    push!(lyaps_Ns,lyaps)
end
pl_pc = plot(ws,critical_curve_adim.(ws);ls=:dash,lc=:gray10,label=L"g_c(w)");
vline!(pl_pc,[(Ns[5]-1)/Ns[5]],lc=:gray10,ls=:dash,label=L"(N-1)/N");
plot!(pl_pc,ws,gs,linez=lyaps_Ns[5],lw=2,label="$(Ns[5])",c=cgrad(:seismic, rev = false),clims=(-1/2,1/2))

#=
for (i,N) in enumerate(Ns)
    plot!(pl,ws,gs,linez=lyaps_Ns[i],lw=2,label="$N",c=cgrad(:seismic, rev = false),clims=(-1/2,1/2))
end
=#