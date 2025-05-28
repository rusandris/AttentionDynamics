using AttentionDynamics
using OrdinaryDiffEq
using Plots
using LaTeXStrings
using LinearAlgebra
using DynamicalSystemsBase
using DelimitedFiles
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/timeseries_phase_portrait.jl")

include("bif_floq_od_method_parameters.jl")
#run necessary mkpaths
#create dir for figs
figs_dir = "../../paper_figs/"
mkpath(figs_dir)


#----------------------------------1D-----------------------------------
#dynamics 

diffeq = (alg = Rodas5(),dt = 1e-3,adaptive=true,reltol = 1e-8, abstol = 1e-10,maxiters=1e15)

N = 1
g = 1.0
w = 0.8 

#generate standard init cond
#u0 = rand(2N)
u0 = generate_initial_condition(N)

w_values = [0.3,0.495]

p = [g,w,N]

function att_din_single_topic!(du,u,p,t)
    n = convert(Int32,p[3])
    du[1] = p[1]*u[1]*(1 - (1 - p[2])*u[2])
    du[2] = u[1] - u[2]
end

#plot eigenvalues of [1,1]
#-(1+gw) \pm \sqrt{(1+gw)^2 - 4g}
#eigenvalue not correct? TODO CHECK 
#1/2(-1 \pm \sqrt{(1 - 4g} #correct
λ12(g,w) = -(1+g*w) + sqrt( Complex((1+g*w)^2 - 4*g))
#x value at fixed point
#x^{*}(w,M,N) = \frac{1}{1 + \frac{wM}{N-1} - \frac{wN}{N-1}}, \quad M \neq 0.
x_fp(w) = #M=1,N=1 

ws = [0:0.01:1.0;]
eigs1 = λ12.(1.0,ws)
plot(ws,real.(eigs1),lw=1,lc=:gray10)
plot!(ws,imag.(eigs1),lw=1,lc=:gray10)

eigs5 = λ12.(5.0,ws)
plot!(ws,real.(eigs5),lw=1,lc=:red)
plot!(ws,imag.(eigs5),lw=1,lc=:red)

eigs50 = λ12.(50.0,ws)
plot!(ws,real.(eigs50),lw=1,lc=:blue)
plot!(ws,imag.(eigs50),lw=1,lc=:blue)

#create system
ds = CoupledODEs(att_din_single_topic!,u0,p;diffeq=diffeq)
#set_parameter!(ds,2,0.8)
#set_parameter!(ds,1,20)

set_parameter!(ds,1,1.0)
set_parameter!(ds,2,0.5)
traj1,t = trajectory(ds,10,u0;Δt=0.001)
plot(traj1[:,1],traj1[:,2],lw=2,lc=:gray10,legendfontsize=14)

set_parameter!(ds,1,5.0)
set_parameter!(ds,2,0.5)
traj5,t = trajectory(ds,10,u0;Δt=0.001)
plot!(traj5[:,1],traj5[:,2],lw=2,lc=:red,legendfontsize=14)

set_parameter!(ds,1,5.0)
set_parameter!(ds,2,0.75)
traj5,t = trajectory(ds,10,u0;Δt=0.001)
plot!(traj5[:,1],traj5[:,2],lw=2,lc=:red,la=0.5,legendfontsize=14)


set_parameter!(ds,1,50.0)
set_parameter!(ds,2,0.5)
traj50,t = trajectory(ds,10,u0;Δt=0.001)
plot!(traj50[:,1],traj50[:,2],lw=2,lc=:blue,legendfontsize=14)


set_parameter!(ds,1,50.0)
set_parameter!(ds,2,0.1)
traj50,t = trajectory(ds,10,u0;Δt=0.001)
plot!(traj50[:,1],traj50[:,2],lw=2,lc=:blue,la=0.5,legendfontsize=14)