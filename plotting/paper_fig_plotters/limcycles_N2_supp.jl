using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using LinearAlgebra
using StatsBase
#using Attractors
using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

x_fp(w, M, N) = (N - 1) / (N - 1 + w * M - w * N) #fixed points

#----------------------------parameters, settings-----------------------
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 2
g = 5.0
w = 0.58
p = [g,w,N]
wc1 = (N-1)/g
wc2 = (N-1)/N
T_short = 12
T_long = 160
Ts = [T_short, T_long]
Ttr = 500
Δt = 1e-3

#create ds
u0 = rand(2N) 
ds = create_ds(u0,p;diffeq=diffeq)

#two w values 
Δw1 = 0.01
Δw2 = 0.005
Δws = [Δw1, Δw2]
ws = [wc1 + Δw1, wc2 - Δw2]
ws_rounded = round.(ws,digits=3)

#loop over w values
plots_limc = []
xlimits = [(Ttr,12),(0,120)]
ylimits = [[(0, 2), (0, 3), (0,2)],[(0, 5), (0, 8) ,(-0.2,5)]]

xlimits_proj = [(-0.2,2.5),(-0.2,4)]
ylimits_proj = [(-0.2,2.5), (-0.2,4)]
lw = 3
titles = [L"w_c = %$(wc1), w = %$(ws_rounded[1])",
    L"w_d = %$(wc2), w = %$(ws_rounded[2])"]

logocolors = Colors.JULIA_LOGO_COLORS
colors = [logocolors.green logocolors.purple]
xticks_ts = [[Ttr, Ttr + 10], [Ttr:50:Ttr+T_long;]]

#make single combined projection plot
#pl_proj = plot()

lw_proj = 3
colors_projs = [:gray10, :gray40]
ticks_projs = [[0,1,2],[0,1,2,3]] 

for (i,w) in enumerate(ws) 

    ylims = ylimits[i]
    w_rounded = round(w; digits=3)
    set_parameter!(ds,2,w)

    #generate all fixed winner-takes-all points (with M = 1)
    x_wta_val = x_fp(w,1,N)
    x_wta = [x_wta_val,x_wta_val,zeros(2N-2)...]
    fixedpoints_wta = [zeros(2N) for i in 1:N]

    for (i,fp) in enumerate(fixedpoints_wta)
        fp[2*i-1:2*i] .= x_wta_val
    end

    #calculate trajectory
    traj,t = trajectory(ds,Ts[i];Ttr=Ttr,Δt=Δt)

    #calculate distances for every fixed point along traj
    dists_fps = []
    for (i,fp) in enumerate(fixedpoints_wta)
        dists = [norm(traj[i] .- fp) for i in 1:length(traj)]
        push!(dists_fps,dists)
    end

    title = titles[i]
    #plot timeseries 
    labels = [L"x_1" L"x_2" L"x_3"]
    x_vars = Matrix(traj[:][:,1:2:2N]) 
    pl_ts = plot(t, x_vars, ylims=ylims[1],lw=lw,
        xformatter=:auto, ylabel=L"x(t)", xlabel=L"t",
        labels=labels, lc=colors, xticks=xticks_ts[i],
        title=title, left_margin=8Plots.mm)

    i == 2 && hline!(pl_ts, [x_wta_val], lc=:gray10, lw=3, ls=:dot, label="")

    #plot the distance to fixed points along the trajectory
    fp_labels = [L"d_1(t)" L"d_2(t)"]
    pl_d = plot(t, stack(dists_fps), ylims=ylims[2],
        xlabel=L"t", ylabel=L"d_i(t)", xticks=xticks_ts[i],
        labels=fp_labels, lw=lw,title=title)

    #plot trajectory projection
    
    #plot trajectory projection
    pl_proj = plot(traj[:, 1], traj[:, 3],
        xlims=xlimits_proj[i], ylims=ylimits_proj[i], 
        ylabel=L"x_2", xlabel=L"x_1",
        xformatter=:auto, labels=labels, lw=lw_proj, la=1.0,
        legend=true, size=(600, 600), lc=colors_projs[1],
        left_margin=-5Plots.mm, aspect_ratio=1,label=L"w=%$(w)",
        title=titles[i],ticks=ticks_projs[i]) #
    
    #plot fixed points 
    
    #coexistence
    x_co_val = x_fp(w, N, N)
    plot!([x_co_val], [x_co_val], st=:scatter, mc=:orange,
        ms=6,markerstrokecolor=:orange,label=L"x_c^{*}")

    #WTA
    plot!([x_wta_val], [0.0], st=:scatter, 
        mc=:white, ms=8, markerstrokecolor=:blue,label=L"x_d^{*}")
    plot!([0.0], [x_wta_val], st=:scatter, 
        mc=:white, ms=8, markerstrokecolor=:blue,label="")

    #=
    pl_3d = plot(traj[:, 1], traj[:, 3], traj[:, 5],  
        ylabel=L"x_2", xlabel=L"x_1", zlabel=L"x_3",
        xformatter=:auto, labels=labels, lw=4, la=1.0, legend=false,
        aspect_ratio=1, lc=:gray40, title=title, proj_type=:ortho,
        cam=(160,30),left_margin=-25Plots.mm,
        tickfontsize=1)
    =#

    push!(plots_limc, pl_ts, pl_d,pl_proj)

end

#=
#add lens
lens!(pl_proj, [0.5, 1.6], [0.5, 1.6], inset=(1, bbox(0.4, 0.1, 0.5, 0.5)))
plot!(pl_proj[2], xformatter=:none, yformatter=:none)
#add single projection plot
push!(plots_limc, pl_proj)
=#

#hline!(plots[4], [x_sw_val], lc=:gray10, lw=2 ,ls=:dash,label="")
pl_limc = plot(plots_limc..., layout=(2, 3), size=(1200, 600), legend_columns=-1,
tickfontsize=10,guidefontsize = 18, bottom_margin=5Plots.mm, legendfontsize=15, dpi=300)

figs_dir = "../../paper_figs/"
savefig(pl_limc, figs_dir * "distance_to_fixedpoints_N2_g5.png")

#plot!(pl_ts,pl_d,layout=(2,1),size=(400,600),
#    guidefontsize=18,legendfontsize=15,tickfontsize=10)

