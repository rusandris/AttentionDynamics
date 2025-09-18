using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using OrdinaryDiffEq
cd(@__DIR__)

competition_strength(w,rc,K,α,N) = w/(1-w) * rc/(α*K*(N-1))
production_rate(g,α) = g*α


#original params
N = 4
rc = rp = 1.0; K = 1.0; α = 5*1e-3

#param values from sup mat 
cs = [0.3,1.0,3.0]

p_orig = [rp,rc,K,cs[1],α]

#rescaled params
w,g = rescaled_params(rp,rc,K,c,α,N)
f1, f2 = rescaling_factors(rp, rc, K, c, α, N)
g_normed = g/(N-1)
p = [g,w,N]

#initialize systems
u0 = rand(2N)
U0 = u0 .* f1
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
diffeq2 = (alg=Rodas5(), dt=1e-4, dtmax=1e-4, adaptive=true, maxiters=1e8)

#diffeq = (alg=Tsit5(), reltol=1e-8, abstol=1e-10)
ds = create_ds(u0,p;diffeq=diffeq)
ds2 = create_ds(u0, p; diffeq=diffeq2)
ds_orig = create_original_ds(U0, p_orig; diffeq=diffeq)

#plot original model time series 
plots_ts = []
for c in cs 
    @show c
    set_parameter!(ds_orig,4,c)
    traj,t = trajectory(ds_orig,4000;Ttr=4000,Δt=1e-2)
    pl = plot(t,Matrix(traj[:,1:2:2N]))
    push!(plots_ts,pl)
end
plot(plots_ts...,layout=(1,3),size=(800,400))