using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)


N = 300
K = 1.0
c = 2.4
r = 10.0
α0 = 0.005
#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

T = 100
Ttr = 1000
Δt = 0.1

#create classic ds
x0 = rand(N) * 0.01
u0 = zeros(2N)
u0[1:2:end] .= x0

ds = create_original_ds(u0, (r/2, r, K, c, α0); diffeq=diffeq)

traj_classic, times_classic = trajectory(ds, T; Ttr=Ttr, Δt=Δt)

trajm_classic = Matrix(traj_classic)[:, 1:2:end]
pl2 = plot(times_classic[1:500], trajm_classic[1:500, :], legend=false)
plot(times_classic[1:500], trajm_classic[1:500, :], legend=false, yscale=:log10, ylimits=(1e-4, 1e1))

#------------------try varied delay system---------------------
#setup varied delays
αs = rand(N) .* 0.01 .+ 1e-3
ds_vd = create_original_varied_delay(u0, (r, r, K, c, αs); diffeq=diffeq)

#test trajectory

traj_vd, times_vd = trajectory(ds_vd, 3T; Ttr=Ttr, Δt=Δt)

trajm_vd = Matrix(traj_vd)[:,1:2:end]
pl_vd = plot(times_vd[1:500], trajm_vd[1:500, :], legend=false,
    title=L"r_p=%$(r),r_c=%$(r),\alpha \in (0.001,0.01)")
pl_vd = plot(times_vd[1:100], trajm_vd[1:100,:], legend=false,yscale=:log10,ylimits=(1e-4,1e1))


#------------------try varied rate system---------------------
rps = rand(6.0:1e-2:20.0, N)
ds_vr = create_original_varied_rate(u0, (rps, r, K, c, α); diffeq=diffeq)

traj_vr, times_vr = trajectory(ds_vr, T; Ttr=Ttr, Δt=Δt)

trajm_vr = Matrix(traj_vr)[:, 1:2:end]
pl_vr = plot(times_vr, trajm_vr, legend=false,
    title=L"r_p \in (6,20) ,r_c=%$(r),\alpha = %$(α)",xlims=(1000,1055))