using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics
using OrdinaryDiffEq
using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using Plots,LaTeXStrings
using LinearAlgebra
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
function get_lifetimes(charts)
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


#set system params 

N = 300

#original params
rp=5.0;rc = 5.0;K = 1.0;α = 5*1e-3;c=1.8
w,g = rescaled_params(rp,rc,K,c,α,N)

p_orig = [rp,rc,K,c,α]

#critical params from the rescaled 
wc1 = (N-1)/(rp/α)
wc2 = (N-1)/N

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
u01 = rand(2N)
u02 = rand(2N)

ds = create_original_ds(u01,p_orig;diffeq=diffeq)

Ttr = 300
T_traj = 500

#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10
c_vals = [2.4:-0.1:1.7;]

op = OrdinalPatterns{N}()
lifetimes = []
for c in c_vals
    @show c
    set_parameter!(ds,4,c)
    traj,ts = trajectory(ds,T_traj;Ttr=Ttr,Δt=10.0)
    xs = traj[:,1:2:2*N]

    #apply OP 
    #sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
    #@show extrema(sts)
    #OP OVERFLOWS FOR LARGE N

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:,end-k+1:end]
    @show charts[1,:] 

    #get lifetimes
    lts = get_lifetimes(charts)
    push!(lifetimes,lts)
end

plots = []
for (i,lts) in enumerate(lifetimes) 

    min,max = extrema(lts)
    bin_edges = range(min-0.1,max+0.1;length=10)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts,bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib = plot(bin_centres, h.weights, lc=:gray20, markershape=:circle,
        xlabel=L"WPE",mc=:gray10, label=L"c=%$(c_vals[i])",dpi=300)

    push!(plots,pl_distrib)

end
plot(plots...,layout=(1,length(c_vals)),xlabel="lifetime",size=(1000,300),
    bottom_margin=5Plots.mm)