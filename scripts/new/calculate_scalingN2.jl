using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 2
g = 5.0
w = 0.58
p = [g, w, N]
wc1 = (N - 1) / g
wc2 = (N - 1) / N

T = 1500
Ttr = 100
Ttr_hopf = 2500
Δt = 1e-2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
ϵ = 0.5

#create ds
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)

#max distance from bifurcation point
Δw_min = 0.001
Δw_max = 0.05
nr_points = 100
#------------------------homoclinic period scaling------------------------------
#params_hc = [0.45:0.001:0.498;]
params_hc = range(wc2 - Δw_max, wc2 - Δw_min; length=nr_points)
periods = period_scaling(ds,params_hc;T=T,Ttr=Ttr,Δt=Δt,u0=u0,ϵ=ϵ)

#--------------------------hopf maximum scaling---------------------------------
#params_hopf = [(N-1)/g:0.001:(N-1)/g+0.1]
params_hopf = range(wc1 + Δw_min, wc1 + Δw_max; length=nr_points)

maximums = amplitude_scaling(ds,params_hopf;T=150,Ttr=Ttr_hopf,Δt=Δt,u0=u0)

#output_dir = "/home/rusandris/Documents/att_din_project/AttentionDynamics/plotting/data/"
savedir = "../../data/scaling_data_N2/"
mkpath(savedir)
writedlm(savedir*"scalingdata_N2_g$g"*".txt",hcat(params_hc,periods,params_hopf,maximums))


