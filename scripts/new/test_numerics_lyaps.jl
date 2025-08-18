using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
cd(@__DIR__)
#inspectdr()
gr()


#check if partial limit cycles really exist
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.3 #w = 0.25

p = [g, w, N]

#search a u0 before transition to limit cycle
T = 50
Ttr = 1000
nr_u0s = 1

u0s = [rand(2N) .* 50 for i in 1:20]

for u0 in u0s[1:2]

    λ1 = lyapunov(ds1, 100; u0=u0, Ttr=Ttr)
    λ2 = lyapunov(ds2, 100; u0=u0, Ttr=Ttr)
    λ2_hp = lyapunov(ds2_hp, 100; u0=u0, Ttr=Ttr)

end
