using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
cd(@__DIR__)

#-----------------------------N=2-------------------------------


N = 2 #N = 5
g = 5.0 #g = 20.0
w = 0.3 #w = 0.25
A = 0.8
p_alle = [g, w, N, A]
p = [g, w, N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)
ds_alle = create_ds_alle(u0, p_alle; diffeq=diffeq)

traj_alle, t = trajectory(ds_alle, 100, u0; Ttr=0, Δt=1e-3)
traj, t = trajectory(ds, 150, u0; Ttr=0, Δt=1e-3)

plot(traj_alle[:, 1], traj_alle[:,2], lc=:gray10, label=L"A = %$A")
plot!(traj[:, 1], traj[:, 3], lc=:red, lw=1, label=L"without Alle")

#------------------------------N=10---------------------------------
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.4 #w = 0.25
A = 1.2

p_alle = [g, w, N, A]
p = [g, w, N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)
ds_alle = create_ds_alle(u0, p_alle; diffeq=diffeq)

traj_alle, t = trajectory(ds_alle, 100, u0; Ttr=0, Δt=1e-3)
traj, t = trajectory(ds, 150, u0; Ttr=0, Δt=1e-3)

plot(traj_alle[:, 1], traj_alle[:, 2], lc=:gray10, label=L"A = %$A")
plot!(traj[:, 1], traj[:, 3], lc=:red, lw=1, label=L"without Alle")