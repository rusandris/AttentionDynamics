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

#create system
ds = CoupledODEs(att_din_single_topic!,u0,p;diffeq=diffeq)
set_parameter!(ds,2,0.8)
set_parameter!(ds,1,20)
traj,t = trajectory(ds,10,u0;Δt=0.001)
plot(t,Matrix(traj[:,:]),lw=2,lc=[:gray10 :gray60],label=[L"x" L"y"],legendfontsize=14)