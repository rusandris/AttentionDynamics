using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles

output_dir = "/home/rusandris/Documents/att_din_project/AttentionDynamics/plotting/data/"

N = 2
g = 5.0
w = 0.2
u0 = rand(2N)
p = [g,w,N]

ds = create_ds(u0,p)

T = 500
Ttr = 100
Δt = 1e-2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
ϵ = 0.5

#homoclinic period scaling
#params_hc = [0.45:0.001:0.498;]
params_hc = range(0.45,0.498;length=50)
periods = period_scaling(ds,params_hc;T=T,Ttr=Ttr,Δt=Δt,state=u0,diffeq=diffeq,ϵ=ϵ)

#hopf maximum scaling
#params_hopf = [(N-1)/g:0.001:(N-1)/g+0.1]
params_hopf = range((N-1)/g,(N-1)/g+0.1;length=50)
maximums = amplitude_scaling(ds,params_hopf;T=50,Ttr=2500,Δt=Δt,state=u0,diffeq=diffeq)

writedlm(output_dir*"scalingdata_N2_g$g"*".txt",hcat(params_hc,periods,params_hopf,maximums))


