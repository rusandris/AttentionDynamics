using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using LinearAlgebra

include("defaults.jl")
include("plot_functions/timeseries_phase_portrait.jl")



#--------------system size and parameters------------------
N = 2
g = 5.0
w = 0.495
p = [g,w,N]
u0 = rand(2*N)

ds = create_ds(u0,p)


Ttr = 150 #transient time
T = 100 
Δt = 0.01

diffeq = (alg=Tsit5(),adaptive=true,dtmax=1e-2,dt=1e-3)

cycle = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = u0,p=p,diffeq = diffeq)

velocity = zeros(2*N)
velocity_magnitudes = zeros(size(cycle)[1])

for i in 1:size(cycle)[1]
	ds.f(velocity,cycle[i,:],p,0) #in-place call to calculate derivatives for every state
	velocity_magnitudes[i] = norm(velocity)
end

velm_ts = plot_timeseries(velocity_magnitudes;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=2,lc=:gray10)
