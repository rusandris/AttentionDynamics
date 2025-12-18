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

function get_lifetimes(trajs,k,Δt_chart)

    N = Int(dimension(ds)/2)
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

    return lifetimes_x,lifetimes_y

end

function plot_distribs(lifetimes,Δt_chart;ps,colors,mode=:probability,nbins = 20,coff=10,scale=:identity,
    xlabel=L"\tau^{(x)}_k \cdot \Delta t",title="lifetime distribution",label_par="")

    pl_distrib = plot()

    for i in 1:length(lifetimes)

        #lifetime distribution
        lts = lifetimes[i][1:end]

        #TODO: REDO BIN_EDGES
        min_lts, max_lts = extrema(lts .* Δt_chart)
        bin_edges = range(min_lts, max_lts; length=nbins)
        bin_width = Float64(bin_edges.step)

        h = fit(Histogram, lts .* Δt_chart, bin_edges)
        h = normalize(h, mode=mode)
        #bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
        #@show bin_edges

        #automatic cutoff: throw out empty bins 
        coff = findfirst(x->iszero(x),h.weights)
        if !(isnothing(coff)) 
            probs = h.weights[1:coff-1]

            #plot together
            plot!(pl_distrib,bin_edges[1:coff-2], probs[1:coff-1], 
                markershape=:circle, scale=scale,st=:scatter,
                mc=colors[i],lc=colors[i],title=title, 
                label=label_par*"$(ps[i])", dpi=300)
        else 
            probs = h.weights

            #plot together
            plot!(pl_distrib,bin_edges[1:end-1], probs, 
                markershape=:circle, scale=scale,st=:scatter,
                mc=colors[i],lc=colors[i],title=title, 
                label=label_par*"$(ps[i])", dpi=300)

        end
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
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)


u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

Ttr = 1000
T_traj = 5000
Δt = 0.5 #sampling rate
#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10

#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

k = 10 #top 10

rs = [6:0.1:15;]
#rs4 = [9.0:1.0:12.0;]

rs4 = [6.0:1.0:9.0;]
set_parameters!(ds,p_orig)
rescaled_params.(rs4, rs4, K, c, α, N)

#precompute trajectories 
trajs = precompute_trajectories(ds,rs4,[1,2],1e-1,15000,1000)

#plot parts of trajectories  
plots_trajs = []

for (i,traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false,title=L"r_p=r_c = %$(rs4[i])")
    pl_y = plot(Matrix(traj[1:500,2:2:16]),legend=false,title=L"r_p=r_c = %$(rs4[i])")

    push!(plots_trajs,pl_x,pl_y)
end
pl_trajs = plot(plots_trajs...,layout=(4,2),size=(600,800),dpi=300)
savefig(pl_trajs,"../../paper_figs/trajs_x_y_change_rp_rc_test.png")

#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [5,10,50]

#=
#test lifetime distrib shape 
#tresholded binning?
lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[2],Δts_chart[1])
lts = lifetimes_x[1]
min, max = extrema(lts .* Δts_chart[1])
bin_edges = range(min, max; length=nbins)
bin_width = Float64(bin_edges.step)

h = fit(Histogram, lts .* Δts_chart[1], nbins=25) #bin_edges)

h = normalize(h, mode=:pdf)
#bin_centres = h.edges[1][1:end-1] .+ 0.5 * h.edges[1].step
plot(h.edges[1][2:end-6], h.weights[2:end-6], st=:scatter, scale=:log10)
plot!([50:1:12200;], lognormal.([50:1:12200;], mean(log10.(lts)), std(log10.(lts))), scale=:log10)
=#

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

nbins = 25
coff_x = 15
coff_y = 13

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=rs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

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
    push!(plots_topk,pl)
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
coff_x = 15
coff_y = 24

#keep k the same, change charting time
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)

    lifetimes_x = lifetimes_x_all_cht[Δt_idx]
    lifetimes_y = lifetimes_y_all_cht[Δt_idx]

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=rs4, colors=colors_x, nbins=nbins, coff=coff_x,
        scale=:log10, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p=r_c = ")

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
    push!(plots_charting, pl)
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

#NUMERICAL METHOD FAILS
#Tsit5 with dt=1e-1 not enough 
#dt = 1e-3 or RODAS5
rs4 = [6.0:1.0:9.0;]
set_parameters!(ds,p_orig)
reinit!(ds)
trajs = precompute_trajectories(ds,rs4,1,1e-1,15000,1000)

#plot parts of trajectories  
plots_trajs = []

for (i,traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false,title=L"r_p = %$(rs4[i]), r_c = %$(rc)")
    pl_y = plot(Matrix(traj[1:12500,2:2:16]),legend=false,title=L"r_p = %$(rs4[i]), r_c = %$(rc)")

    push!(plots_trajs,pl_x,pl_y)
end
pl_trajs = plot(plots_trajs...,layout=(4,2),size=(600,800),dpi=300)
savefig(pl_trajs,"../../paper_figs/trajs_x_y_rp_rc_mismatch_test_change_rp.png")

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end


nbins = 40
coff_x = 4
coff_y = 24

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=rs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_c = %$(rc), r_p = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=rs4,colors=colors_y, nbins=nbins, coff=coff_y, 
        scale=:identity, mode=:pdf,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
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
    push!(plots_topk,pl)
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
savefig(pl_topk,"../../paper_figs/lifetime_distributions_rp_rc_mismatch_change_rp_topk_lowr.png")

#WATCH OUT FOR NUMERICAL METHOD ERRORS
# compare with rp \neq rc (mismatch between production and consumption rate)
#ONLY CHANGE rc
#precompute trajectories 
rs4 = [6.0:1.0:9.0;]
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)
set_parameters!(ds,p_orig)
#set_parameter!(ds, 1,6.0)
trajs = precompute_trajectories(ds,rs4,2,1e-1,15000,1000)

#plot parts of trajectories  
plots_trajs = []

for (i,traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false,title=L"r_p = %$(rp), r_c = %$(rs4[i])")
    pl_y = plot(Matrix(traj[1:5000,2:2:16]),legend=false,title=L"r_p = %$(rp), r_c = %$(rs4[i])")

    push!(plots_trajs,pl_x,pl_y)
end
pl_trajs = plot(plots_trajs...,layout=(4,2),size=(600,800),dpi=300)
savefig(pl_trajs,"../../paper_figs/trajs_x_y_rp_rc_mismatch_test_change_rc.png")

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

nbins = 40
coff_x = 4
coff_y = 24

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=rs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p = %$(rp), r_c = ")

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
    push!(plots_topk,pl)
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

    pl_distrib_x = plot_distribs(lifetimes_x,Δt_chart; 
        ps=rs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"r_p = %$(rp), r_c = ")

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
    push!(plots_charting,pl)
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

#------------params---------------

#gs4 = rs4 ./ α
#αs4 = rp ./ gs4

αs4 = [0.01, 0.0085, 0.0025, 0.00175]
rescaled_params.(rp, rp, K, c, αs4, N)

#loop over memory length
#precompute trajectories 
set_parameters!(ds,p_orig)
trajs = precompute_trajectories(ds,αs4,5,1e-1,15000,1000)

#plot parts of trajectories  
plots_trajs = []

for traj in trajs
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false)
    pl_y = plot(Matrix(traj[1:8500,2:2:16]),legend=false)

    push!(plots_trajs,pl_x,pl_y)
end

pl_trajs = plot(plots_trajs...,layout=(4,2),size=(600,800),dpi=300)
savefig(pl_trajs,"../../paper_figs/trajs_x_y_rp_rc_mismatch_test.png")

#precompute lifetimes 
lifetimes_x_allk = []
lifetimes_y_allk = []

for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])
    push!(lifetimes_x_allk,lifetimes_x)
    push!(lifetimes_y_allk,lifetimes_y)
end

nbins = 25
coff_x = 4
coff_y = 24

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 

    lifetimes_x = lifetimes_x_allk[k_idx]
    lifetimes_y = lifetimes_y_allk[k_idx]

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=αs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"\alpha = ")

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
    push!(plots_topk,pl)
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

    pl_distrib_x = plot_distribs(lifetimes_x,Δt_chart; 
        ps=αs4,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"\alpha = ")

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
    push!(plots_charting,pl)
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

#WATCH OUT FOR NUMERICAL METHOD ERRORS

cs = [1.5:0.3:2.5;]
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)


u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

#precompute trajectories 
set_parameters!(ds,p_orig)
trajs = precompute_trajectories(ds,cs,4,1e-1,15000,1000)

#plot parts of trajectories  
plots_trajs = []

for (i,traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false,title=L"c = %$(cs[i])")
    pl_y = plot(Matrix(traj[1:12500,2:2:16]),legend=false,title=L"c = %$(cs[i])")

    push!(plots_trajs,pl_x,pl_y)
end
pl_trajs = plot(plots_trajs...,layout=(4,2),size=(600,800),dpi=300)
savefig(pl_trajs,"../../paper_figs/trajs_x_y_competition_test.png")


nbins = 20
coff_x = 4
coff_y = 24

lognormal_intervals = [100:1:2000,100:1:22000,100:1:52000]
gaussian_intervals = [1:1:1000,100:1:2700,100:1:2700]

plots_topk = []
for (k_idx,k) in enumerate(ks) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,Δts_chart[1])

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=cs,colors=colors_x,nbins=nbins, coff=coff_x, 
        scale=:log10,mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"c = ")

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
    push!(plots_topk,pl)
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
nbins=20
plots_charting = []
for (Δt_idx, Δt_chart) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x, Δt_chart;
        ps=cs, colors=colors_x, nbins=nbins, coff=coff_x,
        scale=:log10, mode=:pdf,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$(ks[1])" * L",Δt=%$(Δt_chart)",
        label_par=L"c = ")

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
    push!(plots_charting, pl)
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
