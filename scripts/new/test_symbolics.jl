using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics
using OrdinaryDiffEq
using StateTransitionNetworks
using DynamicalSystems
using StatsBase
#test multidim OP after critical point 
#what happens to the symbolic dynamics of w is changed?
#PE,WPE
#average lifetime of in the top k list 

#generate top charts from symbolic timeseries
function get_top_charts(sts,k,encoding)
    
    #decode ordinal patterns: int -> rank pattern
    top_list = zeros(Int64,(length(sts),k))
    for (i,s) in enumerate(sts)
        pattern = decode(encoding,s)
        top_list[i,:] .= pattern[end-k+1:end]
    end
    return top_list

end

#calculate lifetimes by counting consecutive appearances
function get_lifetime(charts)
    tlen,k = size(charts) 

    vars = unique(charts)    

    lifetimes = Int64[]
    #non-optim approach: loop through every var
    for v in vars
        lifetime = 0
        prev_occ = false
        #loop through the charts, every timepoint
        for t in 1:tlen

            #increase lifetime in case of occurrence
            occ = v in charts[t,:]
            occ && (lifetime += 1)

            #save when occurrence streak is over 
            if occ == false && prev_occ == true 
                push!(lifetimes,lifetime)
                lifetime = 0
            end
            prev_occ = occ
        
        end

        #if end of time series ended a streak
        lifetime != 0 && (push!(lifetimes,lifetime))

    end

    return lifetimes

end


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

#codify works as intended on multidim data on 3.6.6 but not on 3.8.1?

sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
perms = [decode(op.encoding,sts[i]) for i in 1:length(sts)]

for w in w_vals
    @show w
    set_parameter!(ds,2,w)

    #trajectories
    traj,t = trajectory(ds,500,u01;Ttr=100,Δt=1e-1)
    xs = traj[:,1:2:2*N]

    sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
    P,Q,ρ = calculate_transition_matrix(sts)

    @show bit_number_measures(ρ)

    @show lyapunov(ds,500;u0=u01,Ttr=10000)
    @show information(Shannon(;base=ℯ),op,xs)
    @show information(Shannon(;base=ℯ),wop,xs)
end
