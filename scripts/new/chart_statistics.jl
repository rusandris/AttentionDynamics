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

    @show N,T

    lifetimes_x = []
    lifetimes_y = []

    #loop over precomputed trajectories 
    #calculate chart lifetimes on instantaneous and cumulated activity  

    for traj in trajs
        xs = traj[1:Δt_chart:T, 1:2:2*N]
        ys = traj[1:Δt_chart:T, 2:2:2*N]

        @show length(xs),length(ys)

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

function plot_distribs(lifetimes,Δt_chart;ps,colors,nbins = 20,coff=10,scale=:identity,xlabel=L"\tau^{(x)}_k \cdot \Delta t",title="lifetime distribution",label_par="")

    pl_distrib = plot()

    for i in 1:length(lifetimes)

        #lifetime distribution
        lts = lifetimes[i][1:end]
        min, max = extrema(lts .* Δt_chart)
        bin_edges = range(min, max; length=nbins)
        bin_width = Float64(bin_edges.step)

        h = fit(Histogram, lts .* Δt_chart, bin_edges)#nbins=50)
        h = normalize(h, mode=:probability)
        bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width

        #plot together
        plot!(pl_distrib,bin_centres[1:coff], h.weights[1:coff], 
            markershape=:circle, scale=scale,
            mc=colors[i],lc=colors[i],title=title, 
            label=label_par*"$(ps[i])", dpi=300)

    end

    plot!(pl_distrib, xlabel=xlabel)

    return pl_distrib
end


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
u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

Ttr = 1000
T_traj = 5000
Δt = 0.5 #sampling rate
#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10

#------------params---------------

#gs4 = rs4 ./ α
#αs4 = rp ./ gs4

αs4 = [0.01, 0.0085, 0.0025, 0.00175]
rescaled_params.(rp, rp, K, c, αs4, N)

#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#loop over memory length
#precompute trajectories 
trajs_αs = precompute_trajectories(ds,αs4,5,0.5,5000,1000)

#plot parts of trajectories  
plots_trajs = []

for traj in trajs
    pl_x = plot(Matrix(traj[1:Δts_chart[1]:500,1:2:2N]),legend=false)
    pl_y = plot(Matrix(traj[1:8500,2:2:4]),legend=false)

    push!(plots_trajs,pl_x,pl_y)
end

#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [5,10,50]

colors_x = palette([:purple, :red], 4)
colors_y = palette([:green, :blue], 4)
plots = []
nbins = 15
coff_x = 7
coff_y = 14
for k in ks 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs_αs,k,1)

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; 
        ps=αs4,colors=colors_x,nbins = nbins, coff=coff_x,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"\alpha = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; 
        ps=αs4,colors=colors_y,nbins = nbins, coff=coff_y,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"\alpha = ")

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)
    push!(plots,pl)
end

plots_charting = []
for (i,Δt_chart) in enumerate(Δts_chart) 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[2],Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x,Δt_chart; 
        ps=rs4,colors=colors_x,nbins = nbins, coff=coff_x,scale=:log10,
        xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"\alpha = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; 
        ps=rs4,colors=colors_y,nbins = nbins, coff=coff_y,scale=:identity,
        xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$(ks[2])" * L",Δt=%$(Δt_chart)",
        label_par=L"\alpha = ")

    @show mean(lifetimes_y[1])
    gaussian = gauss.([100.0:1200;] , mean(lifetimes_y[1])*Δt_chart,std(lifetimes_y[1].*Δt_chart)) 
    plot!(pl_distrib_y,[100.0:1200;],gaussian,lc=:gray10)

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)
    push!(plots_charting,pl)
end

#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

k = 10 #top 10

rs = [6:0.1:15;]
#rs4 = [9.0:1.0:12.0;]

rs4 = [6.0:4.0:18.0;]
set_parameters!(ds,p_orig)
rescaled_params.(rs4, rs4, K, c, α, N)

#precompute trajectories 
trajs = precompute_trajectories(ds,rs4,[1,2],1e-1,5000,1000)

#plot parts of trajectories  
plots_trajs = []

for traj in trajs
    pl_x = plot(Matrix(traj[1:500,1:2:2N]),legend=false)
    pl_y = plot(Matrix(traj[1:500,2:2:2N]),legend=false)

    push!(plots_trajs,pl_x,pl_y)
end

#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [5,10,50]

colors_x = palette([:purple, :red], 4)
colors_y = palette([:green, :blue], 4)
plots_topk = []
for k in ks 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,k,10)

    pl_distrib_x = plot_distribs(lifetimes_x,Δts_chart[1]; ps=rs4,colors=colors_x,nbins = 20,
        coff=10,xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δts_chart[1]; ps=rs4,colors=colors_y,nbins = 20,
        coff=10,xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)
    push!(plots_topk,pl)
end

plots_charting = []
for Δt_chart in Δts_chart 
    lifetimes_x,lifetimes_y =  get_lifetimes(trajs,ks[1],Δt_chart)

    pl_distrib_x = plot_distribs(lifetimes_x,Δt_chart; ps=rs4,colors=colors_x,nbins = 20,
        coff=10,xlabel=L"\tau^{(x)}_k \cdot \Delta t",
        title=L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

    pl_distrib_y = plot_distribs(lifetimes_y,Δt_chart; ps=rs4,colors=colors_y,nbins = 20,
        coff=10,xlabel=L"\tau^{(y)}_k \cdot \Delta t",
        title = L"k=%$k" * L",Δt=%$(Δts_chart[1])",
        label_par=L"r_p=r_c = ")

    pl = plot(pl_distrib_x,pl_distrib_y,layout=(1,2),size=(800,400),
        bottom_margin=5Plots.mm) #yaxis=:log,xaxis=:log)
    push!(plots_charting,pl)
end

plplots



