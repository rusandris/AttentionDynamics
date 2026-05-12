using AttentionDynamics
using DynamicalSystems
using DelimitedFiles
using OrdinaryDiffEq
cd(@__DIR__)

Ns = [4:6;]
w_nr = 10
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2)
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g = gs[30]
Λ = zeros(w_nr,length(Ns))

d0 = 1e-9
T = 5000
Ttr = 10^5

for (n,N) in enumerate(Ns)
    g_N = g * (N-1)
    w_c = (N-1)/g_N
    ws = range(w_c+0.05,(N-1)/N-0.05;length=w_nr)
    u0 = rand(2*N)
    ds = create_ds(u0,[g_N,ws[1],N];diffeq=diffeq)
    states = [u0,u0 .+ d0/sqrt(length(u0))]
    pds = ParallelDynamicalSystem(ds,states)
    @show N
    for (i,w) in enumerate(ws)
        @show w
        @show current_parameters(pds)
        #@show current_states(pds)
        reinit!(pds)
        set_parameter!(pds,2,w)
        Λ[i,n] = lyapunov(pds,T;Ttr=Ttr,Δt=1.0)
        @show Λ[i,n]
    end

end
