using DynamicalSystems,AttentionDynamics,OrdinaryDiffEq,Plots,LaTeXStrings

include("plot_functions/timeseries_phase_portrait.jl")
include("defaults.jl")

p = [6.0,1.0,4]
u0 = rand(8)
ds = create_ds(u0,p)

T = 1.05
Ttr = 0
Δt = 1e-2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=false,maxiters=1e15)
N = 4


ts = evolve_system(ds,T;Ttr = Ttr,Δt = 1e-2,state = u0,p,diffeq=diffeq)

plot_ts = plot_timeseries(ts;T=T,Ttr=Ttr,Δt=Δt,lw=1.5,
xlims = (Ttr,Ttr+30),widen=true)
