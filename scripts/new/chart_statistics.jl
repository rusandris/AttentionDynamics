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


#TODO:
#try exponential/statistically significant binning 
#change charting time (3plot) + write time sampling on axis
#charting time 1/10th of average peak width
#add fits
#show inequality on fig


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


function precompute_trajectories(ds,ps,p_idx,Δt=1e-1,T_traj=5000,Ttr=1000)

    #precomputed trajs
    trajs = []
    for p in ps

        if p_idx isa Vector
            for i in p_idx 
                set_parameter!(ds, i, p)
            end
        else 
            set_parameter!(ds, p_idx, p)
        end

        #solve systems to get trajectory
        @time traj, times = trajectory(ds, T_traj; Ttr=Ttr, Δt=Δt)
        push!(trajs,traj)
    end 
    return trajs
end

function get_lifetimes(ds::DynamicalSystem,ps,p_idxs,Δt,k=10,T_traj=5000,Ttr=1000)

    N = Int(dimension(ds)/2)

    trajs = precompute_trajectories(ds,ps,p_idxs,Δt)

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

    return lifetimes_x,lifetimes_y

end

function get_lifetimes(trajs::Vector{StateSpaceSet}, k, Δt_chart) 

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
 
    N = Int(size(trajs[1])[2]/2)
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
        charts_x = stack([sortperm(xs[i,:]) for i in 1:l])'[:, end-k+1:end]
        charts_y = stack([sortperm(ys[i,:]) for i in 1:l])'[:, end-k+1:end]

        #get lifetimes
        lts_x = get_lifetimes(charts_x)
        lts_y = get_lifetimes(charts_y)
        push!(lifetimes_x, lts_x)
        push!(lifetimes_y, lts_y)
    end

    return lifetimes_x,lifetimes_y

end

function plot_distribs(lifetimes,Δt_chart;ps,colors,mode=:probability,nbins = 50,scale=:log,
    xlabel=L"\tau^{(x)}_k \cdot \Delta t",title="lifetime distribution",label_par="")

    pl_distrib = plot()

    for i in 1:length(lifetimes)

        #lifetime distribution
        lts = lifetimes[i][1:end]

        min_lts, max_lts = extrema(lts .* Δt_chart)

        #apply exponential binning
        b = max_lts ^ (1/nbins) #base 
        bin_edges = min_lts .* b .^ [0:1:nbins;]
        #bin_width = Float64(bin_edges.step)

        h = fit(Histogram, lts .* Δt_chart, bin_edges)
        h = normalize(h, mode=mode)
        #bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
        #@show bin_edges

        #automatic cutoff: throw out empty bins 
        #coff = findfirst(x->iszero(x),h.weights)
        probs = h.weights

        if scale==:log
            X = log10.(bin_edges[1:end-2])
            Y = log10.(probs[1:end-1])
        elseif scale==:lin
            X = bin_edges[1:end-2]
            Y = probs[1:end-1]
        end

        #plot together
        plot!(pl_distrib,X,Y, 
            markershape=:circle,st=:scatter,
            mc=colors[i],lc=colors[i],title=title, 
            label=label_par*"$(ps[i])", dpi=300)

    end

    plot!(pl_distrib, xlabel=xlabel)

    return pl_distrib
end

#calculate lifetimes by counting consecutive appearances
function get_lifetimes(charts;exclude_one=true)
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
        lifetime != 0 && (push!(lifetimes, lifetime))

    end

    return lifetimes

end

#calculate a measure of inequality at each charting time 
#average distance (difference in absolute value) between the first and last position 

function get_inequality(trajs,k,Δt_chart)

    N = Int(dimension(ds)/2)
    T = length(trajs[1]) 

    ineq_avgs = []
    for traj in trajs
        #sampling with charting time 
        xs = traj[1:Δt_chart:T, 1:2:2*N]
        ys = traj[1:Δt_chart:T, 2:2:2*N]

        x_ineqs = Float64[]
        y_ineqs = Float64[]
        for t in 1:length(xs)
            #take top k at each t
            x_topk = sort(xs[t];rev=true)[1:k]
            y_topk = sort(ys[t];rev=true)[1:k]

            #take differences
            x_ineq = x_topk[1] - x_topk[end] 
            y_ineq = y_topk[1] - y_topk[end]

            push!(x_ineqs,x_ineq)
            push!(y_ineqs,y_ineq)
        end
        #take mean and push to vec
        push!(ineq_avgs,(mean(x_ineqs),mean(y_ineqs)))
    end
    return ineq_avgs
end

lognormal(x,μ,σ) = 1/(x*σ*sqrt(2π)) * exp(-(log10(x) - μ)^2/(2*σ^2)) 


function perm_entropy(charts)


end

gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2))

#change 1/alpha memory length
#and see average lifetimes change
αs = [2*1e-4:1e-4:7*1e-3;]




#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

rs4 = [6.0:1.0:9.0;]
Δt = 1e-1
rc = 9.0;
rp = 9.0;
K = 1.0;
c = 2.4;
α = 0.005;

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rp=rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end


#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [5,10,50]

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

#plot lifetimes 

#colors_x = palette([:purple, :red], 4)
#colors_y = palette([:green, :blue], 4)
colors_x = distinguishable_colors(4)
colors_y = distinguishable_colors(4)

nbins = 40

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#
    

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end

xlimits = [(50,1000),(500,10000),(1000,50000)]
for i in 1:3
    xlims!(plots_topk[i][1],xlimits[i])
end 
plot!(plots_topk[3][1],legend=:bottomleft)
plot!(plots_topk[2][1],legend=:bottomleft)
plot!(plots_topk[3][2],legend=:topright)

pl_topk = plot(plots_topk...,layout=(3,1),size=(600,800),dpi=300)

savefig(pl_topk,"../../paper_figs/lifetime_distributions_r_topk_shapes_recolor_lowr.png")


#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx,Δt) in enumerate(Δts_chart) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[2],Δt)
    push!(lifetimes_x_all_cht,lifetimes_x)
    push!(lifetimes_y_all_cht,lifetimes_y)
end

nbins = 40

#keep k the same, change charting time
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_all_cht[Δt_idx]
    lifetimes_y = lifetimes_y_all_cht[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[2])" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=rs4, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p=r_c = ")

    #plot normal distributions 
    for (i, lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts) * Δt_chart, std(lts .* Δt_chart))
        plot!(pl_distrib_y, gaussian_intervals[Δt_idx], gaussian, lc=colors_y[i], label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x, pl_distrib_y, layout=(1, 2), size=(800, 400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end

xlimits = [(500, 5000), (500, 5000), (50, 8000)]
for i in 1:3
    xlims!(plots_charting[i][1], xlimits[i])
end
plot!(plots_charting[3][1], legend=:bottomleft)
plot!(plots_charting[2][1], legend=:bottomleft)
plot!(plots_charting[3][2], legend=:topright)

pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), dpi=300)

savefig(pl_ch, "../../paper_figs/lifetime_distributions_r_charting_time_shapes_recolor_lowr.png")


#-------------------------plot phase diagram---------------------
#check lifetime distributions while keeping rp/alpha constant

#plot critical curve
ws = [0.01:0.01:0.99;]
pl_crit_curve = plot(ws, 1 ./ ws, lc=:gray10, lw=1, ylabel=L"g/(N-1)", xlabel=L"w",
    guidefontsize=15, label=L"w^H",title = L"N=%$(N), K=%$(K), c=%$(c)",legendfontsize=15)

r_test = 9.0

#when rp=rc
w_g_rs4 = rescaled_params.(r_test, r_test, K, c, α, N)
w_rs4 = [w_g_rs4[i][1] for i in 1:length(w_g_rs4)]
g_rs4 = [w_g_rs4[i][2] for i in 1:length(w_g_rs4)]

#plot rp = rc = r params 
plot!(pl_crit_curve, w_rs4, g_rs4 ./ (N - 1),
    mc=:green,lc=:green, lw=1, label=L"r_p=r_c=r \in [6,10,14,18]",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

#when rp \neq rc, only rp changes
w_g_rp = rescaled_params.(rs4, r_test, K, c, α, N)
w_rp = [w_g_rp[i][1] for i in 1:length(w_g_rp)]
g_rp = [w_g_rp[i][2] for i in 1:length(w_g_rp)]

#plot rp = rc = r params 
plot!(pl_crit_curve, w_rp, g_rp ./ (N - 1),
    mc=:purple, lc=:purple, lw=1, label=L"r_p \in [6,10,14,18], r_c=%$(r_test)",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

#when rp \neq rc, only rc changes
w_g_rc = rescaled_params.(r_test, [6,7,8,9], K, c, α, N)
w_rc = [w_g_rc[i][1] for i in 1:length(w_g_rc)]
g_rc = [w_g_rc[i][2] for i in 1:length(w_g_rc)]

#plot rp = rc = r params 
plot!(pl_crit_curve, w_rc, g_rc ./ (N - 1),
    mc=:blue, lc=:blue ,lw=1, label=label=L"r_p=%$(r_test) , r_c \in [6,7,8,9]",ms=4,
    xlims=(0.0,0.5),ylims=(0.,25),markershape=:circle)

savefig(pl_crit_curve,"../../paper_figs/")

# compare with rp \neq rc (mismatch between production and consumption rate)

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rp/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end


nbins = 40

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=rs4,colors=colors_x,nbins=nbins,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins,mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",scale=:lin,
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end

xlims!(plots_topk[1][2],(0,800))
xlims!(plots_topk[2][2],(0,1200))
xlims!(plots_topk[3][2],(0,2800))
plot!(plots_topk[3][1],legend=:bottomleft)
plot!(plots_topk[2][1],legend=:topright)
plot!(plots_topk[1][1],legend=:topright)
plot!(plots_topk[3][2],legend=:topright)

pl_topk = plot(plots_topk...,layout=(3,1),size=(600,800),dpi=300,legendfontsize=5)
subp1 = plot(plots_topk[1][2])
subp2 = plot(plots_topk[2][2])
subp3 = plot(plots_topk[3][2])
pl_topk = plot(subp1,subp2,subp3,layout=(3,1),size=(300,800),dpi=300)
savefig(pl_topk,"../../paper_figs/lifetime_distributions_rp_rc_mismatch_change_rp_topk_lowr.png")

#precompute lifetimes 
lifetimes_x_all_cht = []
lifetimes_y_all_cht = []

for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[2], Δt)
    push!(lifetimes_x_all_cht, lifetimes_x)
    push!(lifetimes_y_all_cht, lifetimes_y)
end

plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_all_cht[Δt_idx]
    lifetimes_y = lifetimes_y_all_cht[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[2])" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=rs4, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p=r_c = ")

    #plot normal distributions 
    for (i, lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts) * Δt_chart, std(lts .* Δt_chart))
        plot!(pl_distrib_y, gaussian_intervals[Δt_idx], gaussian, lc=colors_y[i], label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x, pl_distrib_y, layout=(1, 2), size=(800, 400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end


pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), dpi=300)

#WATCH OUT FOR NUMERICAL METHOD ERRORS
# compare with rp \neq rc (mismatch between production and consumption rate)
#ONLY CHANGE rc
#precompute trajectories 


#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end


#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

nbins = 40

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]


    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p = %$(rp), r_c = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=    
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x[1:2]) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end

xlimits = [(5,10000),(100,10000),(100,50000)]

for i in 1:3
    xlims!(plots_topk[i][1],xlimits[i])
end 

xlims!(plots_topk[1][2],(0,800))
xlims!(plots_topk[3][2],(0,2800))
plot!(plots_topk[3][1],legend=:bottomleft)
plot!(plots_topk[2][1],legend=:topright)
plot!(plots_topk[1][1],legend=:topright)
plot!(plots_topk[3][2],legend=:topright)

pl_topk = plot(plots_topk...,layout=(3,1),size=(600,800),dpi=300,legendfontsize=5)
subp1 = plot(plots_topk[1][2])
subp2 = plot(plots_topk[2][2])
subp3 = plot(plots_topk[3][2])
pl_topk = plot(subp1,subp2,subp3,layout=(3,1),size=(300,800),dpi=300)
savefig(pl_topk,"../../paper_figs/lifetime_distributions_rp_rc_mismatch_change_rc_topk_lowr.png")


#keep k the same, change charting time
plots_charting = []
for (Δt_idx,Δt_chart) in enumerate(Δts_chart) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[1],Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[2])" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")


    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p = %$(rp), r_c = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts)*Δt_chart,std(lts.*Δt_chart)) 
        plot!(pl_distrib_y,gaussian_intervals[Δt_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#

    push!(plots_charting, pl_distrib_x)
end

plot!(plots_charting[3][1],legend=:bottomleft)
plot!(plots_charting[2][1],legend=:bottomleft)
plot!(plots_charting[1][1],legend=:bottomleft)

xlimits = [(10,1000),(10,1000),(10,1000)]

for i in 1:3
    xlims!(plots_charting[i][1],xlimits[i])
end 

pl_ch = plot(plots_charting...,layout=(3,1),size=(600,800),dpi=300)
savefig(pl_ch,"../../paper_figs/lifetime_distributions_rp_rc_mismatch_change_rc_charting time.png")



#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/alpha/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end


#------------params---------------

#gs4 = rs4 ./ α
#αs4 = rp ./ gs4

αs4 = [0.01, 0.0085, 0.0025, 0.00175]

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

nbins = 40

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"\alpha = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#
    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk,pl_distrib_x)
end

xlimits = [(5,1000),(800,10000),(800,50000)]

for i in 1:3
    xlims!(plots_topk[i][1],xlimits[i])
end 

xlims!(plots_topk[1][2],(0,800))
xlims!(plots_topk[3][2],(0,2800))
plot!(plots_topk[3][1],legend=:topright)
plot!(plots_topk[2][1],legend=:topright)
plot!(plots_topk[1][1],legend=:bottomleft)
plot!(plots_topk[3][2],legend=:topright)

pl_topk = plot(plots_topk...,layout=(3,1),size=(600,800),dpi=300)
subp1 = plot(plots_topk[1][2])
subp2 = plot(plots_topk[2][2])
subp3 = plot(plots_topk[3][2])
#pl_topk = plot(subp1,subp2,subp3,layout=(3,1),size=(300,800),dpi=300)
savefig(pl_topk,"../../paper_figs/lifetime_distributions_alpha_memory_topk.png")



lognormal_intervals = [100:1:2000,100:1:2000,100:1:2000]
gaussian_intervals = [1:1:1000,100:1:1200,100:1:1200]

#keep k the same, change charting time
plots_charting = []
for (Δt_idx,Δt_chart) in enumerate(Δts_chart) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[1],Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[2])" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=

    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=αs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"\alpha = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts)*Δt_chart,std(lts.*Δt_chart)) 
        plot!(pl_distrib_y,gaussian_intervals[Δt_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end

plot!(plots_charting[3][1],legend=:bottomleft)
plot!(plots_charting[2][1],legend=:bottomleft)
plot!(plots_charting[1][1],legend=:bottomleft)

xlimits = [(10,1000),(10,1000),(10,1000)]

for i in 1:3
    xlims!(plots_charting[i][1],xlimits[i])
end 

pl_ch = plot(plots_charting...,layout=(3,1),size=(600,800),dpi=300)
savefig(pl_ch,"../../paper_figs/lifetime_distributions_alpha_memory_charting time.png")


#----------------------------change c parameter-------------------------
#increase c -> increase competition 

#READ trajectories
data_dir = "../../data/trajs_precomp_lifetimes/c/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end

cs = [1.5:0.3:2.5;]
nbins = 40

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])

    pl_distrib_x = plot_distribs(lifetimes_x, Δts_chart[1];
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")


    #=
    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=cs,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"c = ")

    #plot normal distributions 
    for (i,lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[k_idx], mean(lts)*Δts_chart[1],std(lts.*Δts_chart[1])) 
        plot!(pl_distrib_y,gaussian_intervals[k_idx],gaussian,lc=colors_y[i],label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[k_idx], mean(log10.(lts.*Δts_chart[1])), std(log10.(lts.*Δts_chart[1])))
        plot!(pl_distrib_x,lognormal_intervals[k_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_topk, pl_distrib_x)
end

xlimits = [(5,1000),(100,10000),(1000,50000)]

for i in 1:3
    xlims!(plots_topk[i][1],xlimits[i])
end 

xlims!(plots_topk[1][2],(0,800))
xlims!(plots_topk[3][2],(0,2000))
plot!(plots_topk[3][1],legend=:bottomleft)
plot!(plots_topk[2][1],legend=:bottomleft)
plot!(plots_topk[1][1],legend=:bottomleft)
plot!(plots_topk[3][2],legend=:topright)

pl_topk = plot(plots_topk...,layout=(3,1),size=(600,800),dpi=300)
subp1 = plot(plots_topk[1][2])
subp2 = plot(plots_topk[2][2])
subp3 = plot(plots_topk[3][2])
pl_topk = plot(subp1,subp2,subp3,layout=(3,1),size=(300,800),dpi=300)
savefig(pl_topk,"../../paper_figs/lifetime_distributions_competition_topk.png")


#keep k the same, change charting time
nbins=40
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t", scale=:log,
        title=L"k=%$(ks[2])" * L",Δt=%$(Δts_chart[1]*Δt)",
        label_par=L"r_c = %$(rc), r_p = ")

    #=
    pl_distrib_y = plot_distribs(lifetimes_y, Δt_chart;
        ps=cs, colors=colors_y, nbins=nbins, coff=coff_y,
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"c = ")

    #plot normal distributions 
    for (i, lts) in enumerate(lifetimes_y)
        gaussian = gauss.(gaussian_intervals[Δt_idx], mean(lts) * Δt_chart, std(lts .* Δt_chart))
        plot!(pl_distrib_y, gaussian_intervals[Δt_idx], gaussian, lc=colors_y[i], label="")
    end

    #=
    #plot lognormal distribution 
    for (i,lts) in enumerate(lifetimes_x) 
        lognorm_vals = lognormal.(lognormal_intervals[Δt_idx], mean(log10.(lts.*Δt_chart)), std(log10.(lts.*Δt_chart)))
        plot!(pl_distrib_x,lognormal_intervals[Δt_idx],lognorm_vals, 
            scale=:log10,lc=colors_x[i],label="")
    end
    =#

    pl = plot(pl_distrib_x, pl_distrib_y, layout=(1, 2), size=(800, 400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)

    =#
    push!(plots_charting, pl_distrib_x)
end

plot!(plots_charting[3][1], legend=:bottomleft)
plot!(plots_charting[2][1], legend=:bottomleft)
plot!(plots_charting[1][1], legend=:bottomleft)

xlimits = [(5, 1000), (5, 1000), (5, 1000)]

for i in 1:3
    xlims!(plots_charting[i][1], xlimits[i])
end

pl_ch = plot(plots_charting..., layout=(3, 1), size=(600, 800), dpi=300)
savefig(pl_ch, "../../paper_figs/lifetime_distributions_competition_charting time.png")
