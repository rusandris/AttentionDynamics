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

function perm_entropy(charts) 
        

end


#set system params 

N = 300

#original params
#rp=3.0;rc = 3.0;K = 1.0;α = 5*1e-3;c=1.8
rc = 9.0;rp = 9.0;K = 1.0;c = 2.4;α = 0.005;
w,g = rescaled_params(rp,rc,K,c,α,N)

p_orig = [rp,rc,K,c,α]

#critical params from the rescaled 
wc1 = (N-1)/(rp/α)
wc2 = (N-1)/N

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01,deepcopy(p_orig);diffeq=diffeq)

Ttr = 500
T_traj = 5500
Δt = 0.5 #sampling rate
#test_traj, = trajectory(ds,50;Ttr=300,Δt=1e-2)

k = 10 #top 10

#change c competition strength 
#and see how difference between 
#the top and bottom of the chart changes
c_vals = [3.4:-0.05:1.4;]

#change 1/alpha memory length
#and see average lifetimes change
αs = [2*1e-4:1e-4:7*1e-3;]
#αs = [4*1e-2:1e-3:8*1e-2;]

rs = [6:0.1:15;]
rs4 = [9.0:1.0:12.0;]

gs4 = rs4 ./ α
αs4 = rp ./ gs4

#check if ratio (w,g) is the same 
w_g_rs4 = rescaled_params.(rs4, rs4, K, c, α, N)
w_g_αs4 = rescaled_params.(rp, rp, K, c, αs4, N)

#--------------------test lyapunov exponent--------------------
n = 10
u0s = [0.1 .* rand(2N) for i in 1:n] #initial conditions 
set_parameters!(ds, p_orig)
λ_rs = []
for r in rs4  
    set_parameter!(ds,1,r)
    set_parameter!(ds,2,r)

    @show r
    p = current_parameters(ds)
    @show p
    @show rescaled_params(p..., N)

    λs = [lyapunov(ds, 1000; u0=u0s[i], Ttr=500) for i in 1:n]
    @show λs
    push!(λ_rs,λs)

end 

set_parameters!(ds,p_orig)
λ_αs = []
for α in αs4
    set_parameter!(ds, 5, α)

    @show α
    p = current_parameters(ds)
    @show p
    @show rescaled_params(p..., N)

    λs = [lyapunov(ds, 1000; u0=u0s[i],Ttr = 500) for i in 1:n]
    @show λs
    push!(λ_αs, λs)

end

pl_lyaps_rs = plot(bar(mean.(λ_rs);yerror=std.(λ_rs)),ylims=(0.0,0.05)) 
pl_lyaps_αs = plot(bar(mean.(λ_αs); yerror=std.(λ_αs), c=:red), ylims=(0.0, 0.05))
plot(pl_lyaps_rs, pl_lyaps_αs,layout=(1,2))
#-------------------test permutation entropy---------------------









#-------------------------plot phase diagram---------------------
#check lifetime distributions while keeping rp/alpha constant

#plot critical curve
ws = [0.01:0.01:0.99;]
pl_crit_curve = plot(ws, 1 ./ ws, lc=:gray10, lw=1, ylabel=L"g/(N-1)", xlabel=L"w",
    guidefontsize=15, label=L"w^H",title = L"N=%$(N), K=%$(K), c=%$(c)",legendfontsize=15)

w_g_αs = rescaled_params.(rp, rc, c, K, αs, N)
w_αs = [w_g_αs[i][1] for i in 1:length(w_g_αs)]
g_αs = [w_g_αs[i][2] for i in 1:length(w_g_αs)]

w_g_rs = rescaled_params.(rs, rs, K, c, α, N)
w_rs = [w_g_rs[i][1] for i in 1:length(w_g_rs)]
g_rs = [w_g_rs[i][2] for i in 1:length(w_g_rs)]

w_g_rs4 = rescaled_params.(rs4, rs4, K, c, α, N)
w_rs4 = [w_g_rs4[i][1] for i in 1:length(w_g_rs4)]
g_rs4 = [w_g_rs4[i][2] for i in 1:length(w_g_rs4)]

colors = palette([:purple, :red], 4)

#plot rp = rc = r curve 
plot!(pl_crit_curve, w_rs, g_rs ./ (N - 1), arrow=:arrow,
    lc=:green, lw=1, label=L"r_p=r_c=r")

plot!(pl_crit_curve,w_rs4, g_rs4 ./ (N - 1), st=:scatter,
    mc=collect(colors), markerstrokewidth=0.0, ms=3, 
    label="",xlims=(0.1,0.5), ylims=(0, 12))

op = OrdinalPatterns{N}()
lifetimes = []
x_vals_diffs = []
plots_ts = []


#loop over production rates
for r in rs4
    @show r
    set_parameter!(ds, 1, r)
    set_parameter!(ds, 2, r)

    #solve systems to get trajectory
    traj, times = trajectory(ds, T_traj; Ttr=Ttr, Δt=Δt)
    xs = traj[:, 1:2:2*N]

    pl = plot(times, Matrix(xs), xlims=(500, 550), legend=false)
    push!(plots_ts, pl)

    #apply OP 
    #sts = [encode(op.encoding,xs[i]) for i in 1:length(xs)]
    #@show extrema(sts)

    #OP pattern-integer conversion OVERFLOWS FOR LARGE N

    #build top charts 
    #charts = get_top_charts(sts,k,op.encoding)
    charts = stack([sortperm(xs[i]) for i in 1:length(xs)])'[:, end-k+1:end]
    @show size(charts)
    #save values?
    #diff between top and bottom of the chart
    x_vals_top_k = [xs[i][charts[i, :]] for i in 1:length(xs)]
    @show size(x_vals_top_k)
    Δxs = [x_vals_top_k[i][end] - x_vals_top_k[i][1] for i in 1:length(x_vals_top_k)]
    push!(x_vals_diffs, Δxs)

    #get lifetimes
    lts = get_lifetimes(charts)
    push!(lifetimes, lts)
end

#plot!(pl, w_αs, g_αs ./ (N - 1), lc=:cyan, lw=3, arrow=:arrow, label="")

plots_distrib = []
for (i, lts) in enumerate(lifetimes)

    min, max = extrema(lts)
    bin_edges = range(min - 0.1, max + 0.1; length=70)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts .* Δt, bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib = plot(bin_centres, h.weights, lc=colors[i], markershape=:circle,
        mc=colors[i],label=L"r=%$(rs4[i])", dpi=300)

    push!(plots_distrib, pl_distrib)

end

plot!(plots_distrib[end],xlabel=L"\tau_k \cdot \Delta t")

pl_distr_r = plot(plots_distrib..., layout=(length(plots_distrib), 1), size=(300, 500),
    bottom_margin=1Plots.mm, xlims=(0, 25), ylims=(0, 0.4), left_margin=8Plots.mm,
    legend=:topright)

pl_r = plot(pl_crit_curve,pl_distr_r,layout=(1,2),size=(800,600),dpi=300)
figs_dir = "../../paper_figs/"
savefig(pl_r,figs_dir*"lifetime_distribution_test_N$(N)_K$(K)_c$(c)_alpha_$(α)_k$(k)"*".png")


#add param annotations to plot
for i in 1:4 
    αs4_round = round(αs4[i];digits=4)
    annotate!(pl_crit_curve,w_rs4[i]-0.01, g_rs4[i] / (N - 1) , text(L"( r=%$(rs4[i]), \alpha = %$(α) )", :gray10,:right,8) )
    annotate!(pl_crit_curve,w_rs4[i] + 0.01, g_rs4[i] / (N - 1) , text(L"( r=%$(rp), \alpha = %$(αs4_round) )", :gray10,:left,8) )
end


lifetimes_alpha = []

#loop over memory length
for α in αs4
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
    push!(lifetimes_alpha,lts)
end

plot(plots_ts...,layout=(length(α_vals),1),size=(500,1600),ylims=(0.0,0.04),left_margin=10Plots.mm)

plots_distrib_alpha = []
for (i,lts) in enumerate(lifetimes_alpha) 

    min, max = extrema(lts)
    bin_edges = range(min - 0.1, max + 0.1; length=70)
    bin_width = Float64(bin_edges.step)
    h = fit(Histogram, lts .* Δt, bin_edges)#nbins=50)
    h = normalize(h, mode=:probability)
    bin_centres = bin_edges[1:end-1] .+ 0.5 * bin_width
    pl_distrib = plot(bin_centres, h.weights, lc=colors[i], markershape=:circle,
        mc=colors[i],label=L"\alpha=%$(αs4[i])", dpi=300)

    push!(plots_distrib_alpha, pl_distrib)

end

plot!(plots_distrib_alpha[end],xlabel=L"\tau_k \cdot \Delta t")

pl_distr_α = plot(plots_distrib_alpha..., layout=(length(plots_distrib_alpha), 1), size=(300, 500),
    bottom_margin=1Plots.mm, xlims=(0, 25), ylims=(0, 0.4), left_margin=8Plots.mm,
    legend=:topright)

pl_α = plot(pl_crit_curve,pl_distr_α,layout=(1,2),size=(800,600),dpi=300)

savefig(pl_α,figs_dir*"lifetime_distribution_test_N$(N)_rp$(rp)_rc$(rp)_K$(K)_c$(c)_k$(k)"*".png")



#loop over competition strength

for c in c_vals
    @show c
    set_parameter!(ds,4,c)

    #solve systems to get trajectory
    traj,times = trajectory(ds,T_traj;Ttr=Ttr,Δt=1.0)
    xs = traj[:,1:2:2*N]

    pl = plot(times,Matrix(xs),xlims=(T_traj,T_traj+100),legend=false)
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

end

plot(plots_ts...,layout=(length(α_vals),1),size=(500,2600),ylims=(0.0,0.04),left_margin=10Plots.mm)
plot(c_vals,mean.(x_vals_diffs),markershape=:circle,ms=1,lc=:gray10,mc=:gray10,
    ylabel=L"\langle x^{(1)} - x^{(k)} \rangle",xlabel=L"c")

#=


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

=#