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
cd(@__DIR__)

#------------------IDEA---------------------
#test multidim OP after critical point 
#what happens to the symbolic dynamics of w is changed?
#PE,WPE
#distribution / average lifetime of in the top k list 
#changing a clearly interpretable parameter: α,c 
#corresponds to moving along a curve in g-w space
#which parameter change leads to driving the system more towards criticality?
#how does properties of the system change when it gets near criticality?

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

#calculate lifetimes by counting consecutive appearances
function get_lifetimes(charts)
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
                push!(lifetimes, lifetime)
                lifetime = 0
            end
            prev_occ = occ

        end

        #if end of time series ended a streak
        lifetime != 0 && (push!(lifetimes, lifetime))

    end

    return lifetimes

end

function perm_entropy(charts)


end

#change 1/alpha memory length
#and see average lifetimes change
αs = [2*1e-4:1e-4:7*1e-3;]


N = 300

#original params
#rp=3.0;rc = 3.0;K = 1.0;α = 5*1e-3;c=1.8
rc = 9.0;
rp = 9.0;
K = 1.0;
c = 2.4;
α = 0.005;
w, g = rescaled_params(rp, rc, K, c, α, N)

p_orig = [rp, rc, K, c, α]

#critical params from the rescaled 
wc1 = (N - 1) / (rp / α)
wc2 = (N - 1) / N

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

Ttr = 1000
T_traj = 500
Δt = 0.5 #sampling rate
#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10


lifetimes_x = []
lifetimes_y = []

x_vals_diffs = []

#loop over memory length
#compute and store trajectories beforehand
trajs = []
for α in αs4
    @show α
    set_parameter!(ds, 5, α)

    #solve systems to get trajectory
    @time traj, times = trajectory(ds, T_traj; Ttr=Ttr, Δt=Δt)
    push!(trajs,traj)
end 

#loop over precomputed trajectories 
#calculate chart lifetimes on instantaneous and cumulated activity  

for traj in trajs
    xs = traj[:, 1:2:2*N]
    ys = traj[:, 2:2:2*N]

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts_x = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:, end-k+1:end]
    charts_y = stack([sortperm(ys[i]) for i in 1:length(ys)])'[:, end-k+1:end]


    #save values?
    #diff between top and bottom of the chart
    #x_vals_top_k = [xs[i][charts[i, :]] for i in 1:length(xs)]
    #@show size(x_vals_top_k)
    #Δxs = [x_vals_top_k[i][end] - x_vals_top_k[i][1] for i in 1:length(x_vals_top_k)]
    #push!(x_vals_diffs, Δxs)

    #get lifetimes
    lts_x = get_lifetimes(charts_x)
    lts_y = get_lifetimes(charts_y)
    push!(lifetimes_x, lts_x)
    push!(lifetimes_y, lts_y)
end

plots_distrib_x = []
plots_distrib_y = []

colors_x = palette([:purple, :red], 4)
colors_y = palette([:green, :blue], 4)

nbins = 30
for i in 1:length(lifetimes_x)

    #plot x lifetime distribution
    lts_x = lifetimes_x[i]
    min, max = extrema(lts_x)
    bin_edges = range(min - 0.1, max + 0.1; length=nbins)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts_x .* Δt, bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib_x = plot(bin_centres, h.weights, markershape=:circle,
        mc=colors_x[i],lc=colors_x[i], label=L"\alpha=%$(αs4[i])", dpi=300)

    push!(plots_distrib_x, pl_distrib_x)

    #plot y lifetime distribution
    lts_y = lifetimes_y[i]
    min, max = extrema(lts_y)
    bin_edges = range(min - 0.1, max + 0.1; length=nbins)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts_y .* Δt, bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib_y = plot(bin_centres, h.weights, markershape=:circle,
        mc=colors_y[i], lc=colors_y[i], label=L"\alpha=%$(αs4[i])", dpi=300)

    push!(plots_distrib_y, pl_distrib_y)
end

plot!(plots_distrib_x[end], xlabel=L"\tau_k \cdot \Delta t")
plot!(plots_distrib_y[end], xlabel=L"\tau_k \cdot \Delta t")

pl_distr = plot(plots_distrib_x..., layout=(length(plots_distrib_x), 1), size=(300, 500),
    bottom_margin=1Plots.mm, left_margin=8Plots.mm,
    legend=:topright)

pl_distr = plot(plots_distrib_y..., layout=(length(plots_distrib_y), 1), size=(300, 500),
    bottom_margin=1Plots.mm, left_margin=8Plots.mm, legend=:topright) #xaxis=:log, yaxis=:log)
    
pl_α = plot(pl_crit_curve, pl_distr_α, layout=(1, 2), size=(800, 600), dpi=300)

savefig(pl_α, figs_dir * "lifetime_distribution_test_N$(N)_rp$(rp)_rc$(rp)_K$(K)_c$(c)_k$(k)" * ".png")
