using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics

#test multidim OP after critical point 
#what happens to the symbolic dynamics of w is changed?
#PE,WPE
#average lifetime of in the top k list 

N = 6
g = 12.0

wc1 = (N-1)/g
wc2 = (N-1)/N
w = 0.5
p = [g,w,N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
u01 = rand(2N)
u02 = rand(2N)

ds = create_ds(u01,p;diffeq=diffeq)

Ttr = 10000
T_traj = 100
T_λ = 2000
n_psos = 5000

w_vals = [0.42,0.45,0.47,0.5]

#test
set_parameter!(ds,2,w_vals[1])
traj,t = trajectory(ds,500;Ttr=1000,Δt=1e-2)
xs = traj[:,1:2:2*N]

op = OrdinalPatterns{N}()
wop = WeightedOrdinalPatterns{N}()
sts = codify(op,xs)

for w in w_vals
    @show w
    set_parameter!(ds,2,w)

    #trajectories
    traj,t = trajectory(ds,500,u01;Ttr=10000,Δt=1e-2)
    xs = traj[:,1:2:2*N]

    @show lyapunov(ds,500;u0=u01,Ttr=10000)
    @show information(op,xs)
    @show information(wop,xs)
end
