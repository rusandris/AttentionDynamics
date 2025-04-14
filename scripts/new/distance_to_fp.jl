using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)

function dist_to_fixedpoint(u,fp)
    return norm(u - fp)
end


#set parameters
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

#N = 5
N = 3
#N = 5 
#g = g_N * (N-1) 
g = 5.0
wc1 = (N-1)/g
wc2 = (N-1)/N
w_test = wc1+ (wc2-wc1)/2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

ds = create_ds(rand(2N),[g,w_test,N];diffeq=diffeq)
traj,t = trajectory(ds,300;Ttr=1000)


#dominance fixed point: doesn't matter which one we choose 
fp_dom = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=1
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)

fp_dom[[1,N+1]] .= [s,s] #set both x1 and y1 to s

dists_dom = [dist_to_fixedpoint(traj[i],fp_dom) for i in 1:length(traj)]

ws = [0.3:0.001:0.9;]

for (i,w) in enumerate(ws)
    dists_dom = [dist_to_fixedpoint(traj[i],fp_dom) for i in 1:length(traj)]


end