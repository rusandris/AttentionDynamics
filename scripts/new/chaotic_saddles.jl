using Revise
using BifurcationKit
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using Plots
using DynamicalSystems


#-----------IDEA---------------
#plot chaotic saddle set 
#plot transient part of trajectories with positive lyap exp

N = 10
g = 45.0

wc1 = (N-1)/g
wc2 = (N-1)/N
w = 0.3
p = [g,w,10]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
u0 = rand(2N)
ds = create_ds(u0,p;diffeq=diffeq)
set_parameter!(ds,2,0.4)
traj,t = trajectory(ds,500;Ttr=0,Δt=1e-2)
lyapunov(ds,2000;u0=u0,Ttr=1000)