export get_lifetimes,get_lifetimes_integ, get_distrib, integrate_activity
export chart_stats, consecutive_occurrences
export autohist, parabolic, linear #plot_distrib,plot_distrib!


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

function get_lifetimes(trajs, k, Δt_chart)

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

function get_lifetimes(traj::Matrix{Float64}, k, Δt_chart)

    T,d = size(traj)
    N = Int(d / 2)


    xs = @view traj[1:Δt_chart:T, 1:2:2*N]
    ys = @view traj[1:Δt_chart:T, 2:2:2*N]
    l = size(xs)[1]
    #@show length(xs),length(ys)

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts_x = stack([sortperm(xs[i, :]) for i in 1:l])'[:, end-k+1:end]
    charts_y = stack([sortperm(ys[i, :]) for i in 1:l])'[:, end-k+1:end]

    #get lifetimes
    lts_x = get_lifetimes(charts_x)
    lts_y = get_lifetimes(charts_y)

    return lts_x, lts_y

end


#TODO: implement statistically significant binning

function get_distrib(τs; binning=:auto, mode=:none, nbins=50, scale=:log)
    
    #decide bin edges automatically
    #number of bins is given
    min_lts, max_lts = extrema(τs)
    if binning ==:exp 
        #apply exponential binning
        b = max_lts^(1 / nbins) #base 
        bin_edges = min_lts .* b .^ [0:1:nbins;]
    elseif binning ==:lin
        bin_edges = range(min_lts, max_lts; length=nbins)
    end
    
    get_distrib(τs, bin_edges; mode=:none, scale=scale)

end

function get_distrib(τs, bin_edges; mode=:none, scale=:log)

    h = fit(Histogram, τs, bin_edges)

    h = normalize(h, mode=mode)
    #bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    #@show bin_edges

    Y = h.weights
    X = h.edges[1][1:end-1]

    any(iszero.(Y)) && @warn "There are bins with zero probability!"
    @show length(unique(Y))
    @show Y
    if scale == :log
        return log10.(X), log10.(Y)
    end

    return X, Y
    
end

#determine exponential bin edges so that each contains at least m points
function autohist(y;m=10,nbins=50,b=10,mode=:none) 
    ymin,ymax = extrema(y) #interval to cover

    n = nbins #start from here
    bins = b .^ range(log(b, ymin), log(b, ymax), length=nbins)
    h = fit(Histogram, y, bins)
    n_min = 0
    while true
        h = fit(Histogram, y, bins)
        n_min = minimum(h.weights)
        n_min >= m && break

        n = n - 1 
        bins = b .^ range(log(b, ymin), log(b, ymax), length=n)
    end
    h = normalize(h, mode=mode)
    return h.edges[1],h.weights
end 

#=
function plot_distrib(τs, Δt_chart; binning=:auto, color=:gray10, st=:scatter, mode=:none, ms=2, nbins=50, scale=:log,
    xlabel=L"\log(\tau^{(x)}_k \cdot \Delta t)", ylabel=L"\log(count)", title="lifetime distribution", label_par="")

    X,Y = get_distrib(τs, Δt_chart; binning=binning, color=color, mode=mode, nbins=nbins, scale=scale)

    #plot together
    plot(X, Y, markershape=:circle, ms=ms, 
        ylabel=ylabel,xlabel=xlabel,
        mc=color, lc=color, title=title,
        label=label_par, dpi=300)

end


function plot_distrib!(pl,τs, Δt_chart; binning=:auto, color=gray10,st=:scatter, ms=2 ,mode=:none, nbins=50, scale=:log,
    xlabel=L"\log(\tau^{(x)}_k \cdot \Delta t)", ylabel=L"\log(count)", title="lifetime distribution", label_par="")

    X, Y = get_distrib(τs, Δt_chart; binning=binning, color=color, mode=mode, nbins=nbins, scale=scale)

    #plot together
    plot!(pl, X, Y, markershape=:circle, ms=ms,
        ylabel=ylabel,xlabel=xlabel,st=st,
        mc=color,markerstrokecolor=color, lc=color, title=title,
        label=label_par, dpi=300)

end
=#

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

function chart_stats(traj, k, Δt_chart;threshold=1e-3)

    T,d = size(traj) 

    #sampling with charting time 
    xs = traj[1:Δt_chart:T, 1:d]

    x_ineqs = Float64[]
    nr_inside_thresh = Int64[]
    x_topk_vals = Vector{Float64}[]
    x_topk_idxs = Vector{Int64}[]
    for t in 1:size(xs)[1]
        #build chart top k at each t
        xt = xs[t, :]
        chart = reverse(sortperm(xt)[end-k+1:end])
        x_topk = sort(xt; rev=true)[1:k]

        #count how many inside threshold
        x_topk_ratios = x_topk ./ x_topk[1]
        inside_thresh = count(x -> x > threshold,x_topk_ratios)

        #take relative differences 
        x_ineq = (x_topk[1] - x_topk[end]) / x_topk[end]

        push!(x_ineqs, x_ineq)
        push!(x_topk_vals,x_topk)
        push!(x_topk_idxs, chart)
        push!(nr_inside_thresh, inside_thresh)
    end
    #take mean and push to vec
    return x_topk_idxs, x_topk_vals, nr_inside_thresh, extrema(x_ineqs), mean(x_ineqs)
end

function consecutive_occurrences(charts)
    idxs_active = unique(stack(charts))
    occs_idxs = []
    for i in idxs_active
        occs = zeros(length(charts))
        for (c,chart) in enumerate(charts)
            rank = indexin(i,chart)[1]
            if isnothing(rank)
                occs[c] = NaN
            else 
                occs[c] = rank
            end
        end
        push!(occs_idxs,occs)
    end
    return occs_idxs
end


lognormal(x, μ, σ) = 1 / (x * σ * sqrt(2π)) * exp(-(log10(x) - μ)^2 / (2 * σ^2))
parabolic(L, a, b) = 10^(-a*log10(L)-b*log10.(L)^2)
linear(L, a, b) = 10^(a * log10(L) + b )

function perm_entropy(charts)


end

function get_lifetimes_integ(Is, k)

    charts = stack([sortperm(Is[i]) for i in 1:length(Is)])'[:, end-k+1:end]

    lts = get_lifetimes(charts)

end

function integrate_activity(traj, ts, n::Int; method=Integrals.TrapezoidalRule())
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