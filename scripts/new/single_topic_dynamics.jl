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

#plot eigenvalues
#-(1+gw) \pm \sqrt{(1+gw)^2 - 4g}
#eigenvalue not correct? TODO CHECK
#eigenvalue of N>1 identical initial conditions
#λ12(g,w) = -(1+g*w) + sqrt( Complex((1+g*w)^2 - 4*g))
#x value at fixed point

#eigenvale for [x^*,x^*] N=1 case
#1/2(-1 \pm \sqrt{(1 - 4g} #correct
λ12(g) = 1/2*(-1 + sqrt( Complex(1 - 4g) ))

#x^{*}(w,M,N) = \frac{1}{1 + \frac{wM}{N-1} - \frac{wN}{N-1}}, \quad M \neq 0.
x_fp(w) = 1/(1-w) #M=1,N=1 

ws = [0:0.01:1.0;]
gs = [2:50;]
eigs = λ12.(gs)
plot(gs,real.(eigs),lw=1,lc=:gray10,xlabel=L"g",guidefontsize=20,label="real")
plot!(gs,imag.(eigs),lw=1,lc=:red,label="imag")

#create system
ds = CoupledODEs(att_din_single_topic!,u0,p;diffeq=diffeq)
#set_parameter!(ds,2,0.8)
#set_parameter!(ds,1,20)

#plot trajectories in 2d
#plot time series?
#show time scale separation role of g 
#fixed point value depends on w 
#fixed point eigenvalue only depends on g 

p1 = [1.0,0.5,1.0]
p100 = [100.0,0.5,1.0]


set_parameters!(ds,p1)
traj1,t = trajectory(ds,10,u0;Δt=0.001)

#plot time series
pl_ts = plot(t,traj1[:,1],lw=2,lc=:gray10,la=1,legendfontsize=14,
    xlabel=L"t",ylabel=L"x",guidefontsize=15,label=L"g=1")

set_parameters!(ds,p100)
traj100,t = trajectory(ds,10,u0;Δt=0.001)

plot!(pl_ts,t,traj100[:,1],lw=2,lc=:blue,la=1,legendfontsize=14,label=L"g=100",ylims=[0,10])

#---------------phase space velocity plots-----------------
#plot under time series?
velocity = zeros(4)
velocity_magnitudes1 = zeros(length(traj1[:,1]))
velocity_magnitudes100 = zeros(length(traj100[:,1]))

for i in 1:size(traj1)[1]
    #in-place call to calculate derivatives for every state
    att_din_single_topic!(velocity,traj1[i],p1,0) 
    velocity_magnitudes1[i] = norm(velocity)
    att_din_single_topic!(velocity,traj100[i],p100,0)
    velocity_magnitudes100[i] = norm(velocity)
end

pl_vs = plot(t,velocity_magnitudes1,lw=2,lc=:gray10,la=1,legendfontsize=14,
    xlabel=L"t",ylabel=L"|v|",guidefontsize=15,label=L"g=1",ylims=(0,10))

plot!(pl_vs,t,velocity_magnitudes100,lw=2,lc=:blue,la=1,legendfontsize=14,
    xlabel=L"t",ylabel=L"|v|",guidefontsize=15,label=L"g=100")

pl2 = plot(pl_ts,pl_vs,layout=(2,1),size=(600,800))

#plot phase portrait
pl_pp = plot(traj1[:,1],traj1[:,2],lw=2,lc=:gray10,legendfontsize=14,label=L"g=1",
    xlabel=L"x",ylabel=L"y",guidefontsize=15)
plot!(pl_pp,traj100[:,1],traj100[:,2],lw=2,lc=:blue,legendfontsize=14,label=L"g=100")

pl = plot(pl2,pl_pp,layout=(1,2),size=(1200,500),
    left_margin=10Plots.mm,bottom_margin=10Plots.mm,
    dpi=300)
savefig(pl,"single_topic_dynamics_test.png")

set_parameter!(ds,1,100.0)
set_parameter!(ds,2,0.9)
traj_test,t = trajectory(ds,10,u0;Δt=0.001)
plot!(traj_test[:,1],traj_test[:,2],lw=2,lc=:blue,la=0.4,legendfontsize=14,label=L"g=100")
plot!(t,traj_test[:,1],lw=2,lc=:blue,la=0.4,legendfontsize=14,label=L"g=100")

#quick anim 
anim = @animate for i in 1:50:length(traj_test)
    plot(traj_test[:,1],traj_test[:,2],lw=2,lc=:blue,
        la=0.4,legendfontsize=14,label=L"g=100",
        xlabel=L"x",ylabel=L"y",guidefontsize=15)
    plot!([traj_test[i,1]],[traj_test[i,2]],st=:scatter,ms=5,mc=:gray10,label="")
end
gif(anim, "anim_fps15.gif", fps = 10)


