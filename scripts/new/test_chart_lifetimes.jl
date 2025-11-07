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
u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01,p_orig;diffeq=diffeq)

Ttr = 500
T_traj = 1500
Δt = 1.0 #sampling rate
#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10

#change c competition strength 
#and see how difference between 
#the top and bottom of the chart changes
c_vals = [2.4:-0.1:1.7;]

#change 1/alpha memory length
#and see average lifetimes change
α_vals = [2*1e-4:1e-4:7*1e-3;]

op = OrdinalPatterns{N}()
lifetimes = []
x_vals_diffs = []
plots_ts = []

#=
#loop over competition strength
for c in c_vals
    @show c
    set_parameter!(ds,4,c)

    #solve systems to get trajectory
    traj,times = trajectory(ds,T_traj;Ttr=Ttr,Δt=1.0)
    xs = traj[:,1:2:2*N]

    pl = plot(times,Matrix(xs),xlims=(450,500),legend=false)
    push!(plots_ts,pl)

    #apply OP 
    #sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
    #@show extrema(sts)

    #OP pattern-integer conversion OVERFLOWS FOR LARGE N

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:,end-k+1:end]
    @show size(charts)
    #save values?
    #diff between top and bottom of the chart
    x_vals_top_k =  [ xs[i][charts[i,:]] for i in 1:length(xs)]
    @show size(x_vals_top_k)
    Δxs = [x_vals_top_k[i][end] -  x_vals_top_k[i][1] for i in 1:length(x_vals_top_k)]
    push!(x_vals_diffs,Δxs)

    #get lifetimes
    lts = get_lifetimes(charts)
    push!(lifetimes,lts)
end
=#

#loop over memory length
for α in α_vals
    @show α
    set_parameter!(ds,5,α)

    #solve systems to get trajectory
    traj,times = trajectory(ds,T_traj;Ttr=Ttr,Δt=1.0)
    xs = traj[:,1:2:2*N]

    pl = plot(times,Matrix(xs),xlims=(500,550),legend=false)
    push!(plots_ts,pl)

    #apply OP 
    #sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
    #@show extrema(sts)

    #OP pattern-integer conversion OVERFLOWS FOR LARGE N

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:,end-k+1:end]
    @show size(charts)
    #save values?
    #diff between top and bottom of the chart
    x_vals_top_k =  [ xs[i][charts[i,:]] for i in 1:length(xs)]
    @show size(x_vals_top_k)
    Δxs = [x_vals_top_k[i][end] -  x_vals_top_k[i][1] for i in 1:length(x_vals_top_k)]
    push!(x_vals_diffs,Δxs)

    #get lifetimes
    lts = get_lifetimes(charts)
    push!(lifetimes,lts)
end

plot(plots_ts...,layout=(length(α_vals),1),size=(500,1600),ylims=(0.0,0.04),left_margin=10Plots.mm)

plots_distrib = []
for (i,lts) in enumerate(lifetimes) 

    min,max = extrema(lts)
    bin_edges = range(min-0.1,max+0.1;length=5)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts,bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib = plot(bin_centres, h.weights, lc=:gray20, markershape=:circle,
        mc=:gray10, label=L"α=%$(α_vals[i])",dpi=300)

    push!(plots_distrib,pl_distrib)

end

#plot(plots_ts...,layout=(length(c_vals),1),size=(500,1600),ylims=(0.0,0.04),left_margin=10Plots.mm)

plot!(plots_distrib[end],xlabel=L"\langle \tau_k \rangle")
plot(plots_distrib...,layout=(length(plots_distrib),1),size=(300,1500),
    bottom_margin=-3Plots.mm,xlims=(0,25),ylims=(0,1),left_margin=8Plots.mm,
    legend=:outertopright)



#plot mean of lifetimes 
pl = plot(α_vals,mean.(lifetimes) .* Δt,ylabel=L"\langle \tau_k \rangle \cdot \Delta t",
    xlabel=L"\alpha",lc=:gray10,mc=:gray10,ms=3,markershape=:circle,
    guidefontsize=15,xaxis=:log,yaxis=:log,dpi=300,
    label = "N=300,Δt=1,k=10",title="Average lifetime inside top k")

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"average_lifetime_test_N$(N)_rp$(rp)_rc$(rp)_K$(K)_c$(c)_k$(k)"*".png")