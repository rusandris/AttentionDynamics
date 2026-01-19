using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics
using OrdinaryDiffEq
using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using Integrals

#generate top charts from symbolic timeseries
function get_top_charts(sts, k, encoding)

    #decode ordinal patterns: int -> rank pattern
    top_list = zeros(Int64, (length(sts), k))
    for (i, s) in enumerate(sts)
        pattern = decode(encoding, s)
        top_list[i, :] .= pattern[end-k+1:end]
    end
    return top_list

end

function get_lifetimes(ds::DynamicalSystem, ps, p_idxs, Δt, k=10, T_traj=5000, Ttr=1000)

    N = Int(dimension(ds) / 2)

    trajs = precompute_trajectories(ds, ps, p_idxs, Δt)

    lifetimes_x = []
    lifetimes_y = []

    #loop over precomputed trajectories 
    #calculate chart lifetimes on instantaneous and cumulated activity  

    for traj in trajs
        xs = traj[:, 1:2:2*N]
        ys = traj[:, 2:2:2*N]

        #build top charts 
        #charts = get_top_charts(sts,k,op.encoding)
        charts_x = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:, end-k+1:end]
        charts_y = stack([sortperm(ys[i]) for i in 1:length(ys)])'[:, end-k+1:end]

        #get lifetimes
        lts_x = get_lifetimes(charts_x)
        lts_y = get_lifetimes(charts_y)
        push!(lifetimes_x, lts_x)
        push!(lifetimes_y, lts_y)
    end

    return lifetimes_x, lifetimes_y

end

function get_lifetimes(trajs::Vector{<:StateSpaceSet}, k, Δt_chart)

    N = Int(dimension(ds) / 2)
    T = length(trajs[1])

    #@show N,T

    lifetimes_x = []
    lifetimes_y = []

    #loop over precomputed trajectories 
    #calculate chart lifetimes on instantaneous and cumulated activity  

    for traj in trajs
        xs = traj[1:Δt_chart:T, 1:2:2*N]
        ys = traj[1:Δt_chart:T, 2:2:2*N]

        #@show length(xs),length(ys)

        #build top charts 
        #charts = get_top_charts(sts,k,op.encoding)
        charts_x = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:, end-k+1:end]
        charts_y = stack([sortperm(ys[i]) for i in 1:length(ys)])'[:, end-k+1:end]

        #get lifetimes
        lts_x = get_lifetimes(charts_x)
        lts_y = get_lifetimes(charts_y)
        push!(lifetimes_x, lts_x)
        push!(lifetimes_y, lts_y)
    end

    return lifetimes_x, lifetimes_y

end

function get_lifetimes(trajs::Vector{Matrix{Float64}}, k, Δt_chart)

    N = Int(size(trajs[1])[2] / 2)
    T = size(trajs[1])[1]

    #@show N,T

    lifetimes_x = []
    lifetimes_y = []

    #loop over precomputed trajectories 
    #calculate chart lifetimes on instantaneous and cumulated activity  

    for traj in trajs
        xs = traj[1:Δt_chart:T, 1:2:2*N]
        ys = traj[1:Δt_chart:T, 2:2:2*N]
        l = size(xs)[1]
        #@show length(xs),length(ys)

        #build top charts 
        #charts = get_top_charts(sts,k,op.encoding)
        charts_x = stack([sortperm(xs[i, :]) for i in 1:l])'[:, end-k+1:end]
        charts_y = stack([sortperm(ys[i, :]) for i in 1:l])'[:, end-k+1:end]

        #get lifetimes
        lts_x = get_lifetimes(charts_x)
        lts_y = get_lifetimes(charts_y)
        push!(lifetimes_x, lts_x)
        push!(lifetimes_y, lts_y)
    end

    return lifetimes_x, lifetimes_y

end

function plot_distribs(lifetimes, Δt_chart; ps, colors, mode=:none, nbins=50, scale=:log,
    xlabel=L"\log(\tau^{(x)}_k \cdot \Delta t)", title="lifetime distribution", label_par="")

    pl_distrib = plot()

    for i in 1:length(lifetimes)

        #lifetime distribution
        lts = lifetimes[i][1:end]

        min_lts, max_lts = extrema(lts .* Δt_chart)

        #apply exponential binning
        b = max_lts^(1 / nbins) #base 
        bin_edges = min_lts .* b .^ [0:1:nbins;]
        #bin_width = Float64(bin_edges.step)

        h = fit(Histogram, lts .* Δt_chart, bin_edges)
        h = normalize(h, mode=mode)
        #bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
        #@show bin_edges

        #automatic cutoff: throw out empty bins 
        #coff = findfirst(x->iszero(x),h.weights)
        probs = h.weights

        if scale == :log
            X = log10.(bin_edges[1:end-2])
            Y = log10.(probs[1:end-1])
        elseif scale == :lin
            X = bin_edges[1:end-2]
            Y = probs[1:end-1]
        end

        #plot together
        plot!(pl_distrib, X, Y,
            markershape=:circle, st=:scatter,ylabel="log(count)",
            mc=colors[i], lc=colors[i], title=title,
            label=label_par * "$(ps[i])", dpi=300)

    end

    plot!(pl_distrib, xlabel=xlabel)

    return pl_distrib
end

#calculate lifetimes by counting consecutive appearances
function get_lifetimes(charts; exclude_one=true)
    tlen, k = size(charts)

    vars = unique(charts)

    lifetimes = Int64[]
    #non-optim approach: loop through every var
    for v in vars
        lifetime = 0
        prev_occ = false
        #loop through the charts, every timepoint
        for t in 1:tlen

            #increase lifetime in case of occurrence
            occ = v in charts[t, :]
            occ && (lifetime += 1)

            #save when occurrence streak is over 
            if occ == false && prev_occ == true

                if exclude_one == true
                    lifetime > 1 && (push!(lifetimes, lifetime))
                else
                    push!(lifetimes, lifetime)
                end
                lifetime = 0
            end
            prev_occ = occ

        end

        #if end of time series ended a streak
        if exclude_one == true
            lifetime > 1 && (push!(lifetimes, lifetime))
        else 
            lifetime != 0 && (push!(lifetimes, lifetime))
        end
        
    end

    return lifetimes

end

#calculate a measure of inequality at each charting time 
#average distance (difference in absolute value) between the first and last position 

function get_inequality(trajs, k, Δt_chart)

    T,d = size(trajs[1]) 
    N = Int(d / 2)

    ineq_avgs = []
    for traj in trajs
        #sampling with charting time 
        xs = traj[1:Δt_chart:T, 1:2:2*N]
        ys = traj[1:Δt_chart:T, 2:2:2*N]

        x_ineqs = Float64[]
        y_ineqs = Float64[]
        for t in 1:T
            #take top k at each t
            x_topk = sort(xs[t,:]; rev=true)[1:k]
            y_topk = sort(ys[t,:]; rev=true)[1:k]

            #take differences
            x_ineq = x_topk[1] - x_topk[end]
            y_ineq = y_topk[1] - y_topk[end]

            push!(x_ineqs, x_ineq)
            push!(y_ineqs, y_ineq)
        end
        #take mean and push to vec
        push!(ineq_avgs, (mean(x_ineqs), mean(y_ineqs)))
    end
    return ineq_avgs
end

lognormal(x, μ, σ) = 1 / (x * σ * sqrt(2π)) * exp(-(log10(x) - μ)^2 / (2 * σ^2))


function perm_entropy(charts)


end


#TODO: 
#implement chart building from integrated variables

function get_lifetimes_integ(Is, k)

    charts = stack([sortperm(Is[i]) for i in 1:length(Is)])'[:, end-k+1:end]

    lts = get_lifetimes(charts)

end

function integrate_activity(traj, ts, n::Int; method=TrapezoidalRule())
    #T = length(traj)
    #l = length(traj[1])
    T,l = size(traj)

    traj_integ = Vector{Float64}[]
    time_segments = []
    start_idxs = [1:n:T;][1:end-1]
    for i in start_idxs
        traj_segment = traj[i:i+n,:]
        t = ts[i:i+n]
        int_vals = integrate_activity(traj_segment, t; method=method)
        push!(traj_integ,int_vals)
        push!(time_segments, (t[1],t[end]))
    end 
    return traj_integ,time_segments
end

function integrate_activity(traj, ts::Vector{Float64}; method=TrapezoidalRule())
    t,nr_vars = size(traj) #nr of variables
    L_ints = zeros(nr_vars)
    for i in 1:nr_vars 
        y = traj[:,i]
        problem = SampledIntegralProblem(y, ts)
        Li_int = solve(problem, method)
        L_ints[i] = Li_int.u
    end
    return L_ints
end