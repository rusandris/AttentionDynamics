using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles

#--------------system size and parameters------------------
N = 4
g = 5.0
w = 0.2
p = [g,w,N]

#----------------------select initial condition----------------------
u0s = readdlm("N4_LIMCYCLES_IC.txt")
ic_index = 5
u0 = u0s[:,ic_index]


d0 = 1e-9 #initial traj distance
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15) 

#create parallel integrator
pinteg = create_pinteg(u0,p,d0;diffeq=diffeq)


Ttr = 15000 #transient time
T = 1000 #comp time

#parameter ranges
gvalues = range(3,10,10)
wvalues = range(0.3,1.0,10)

#phase diagram calculation
lyapunov_diagram(pinteg,gvalues,wvalues;T=T,Ttr=Ttr,Δt = 1,save=true,realtime_plotting=false)
