using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using LinearAlgebra
using StatsBase
#using Attractors
using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)


diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#diffeq = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15)
N = 3
g = 5.0
w = 0.58
p = [g,w,N]
wc1 = (N-1)/g
wc2 = (N-1)/N
T_short = 12
T_long = 120
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

#loop over w values
plots = []
xlimits = [(Ttr,12),(0,120)]
ylimits = [[(0, 3), (0, 3), (0,5)],[(0, 6), (0, 8) ,(0,10)]]

xlimits_proj = [(0,3),(0,6)]
ylimits_proj = [(0,3), (0,6)]
lw = 3
titles = [L"w = w_c + %$Δw1", L"w = w_d - %$Δw2"]

for (i,w) in enumerate(ws) 

    ylims = ylimits[i]
    w_rounded = round(w; digits=3)
    set_parameter!(ds,2,w)

    #generate all fixed single-winner points (with M = 1)
    x_fp(w,M,N) = (N-1)/(N-1 + w*M - w*N)
    x_sw_val = x_fp(w,1,N)
    x_sw = [x_sw_val,x_sw_val,zeros(2N-2)...]
    fixedpoints_sw = [zeros(2N) for i  in 1:N]

    for (i,fp) in enumerate(fixedpoints_sw)
        fp[2*i-1:2*i] .= x_sw_val
    end

    #calculate trajectory
    traj,t = trajectory(ds,Ts[i];Ttr=Ttr,Δt=Δt)

    #calculate distances for every fixed point along traj
    dists_fps = []
    for (i,fp) in enumerate(fixedpoints_sw)
        dists = [norm(traj[i] .- fp) for i in 1:length(traj)]
        push!(dists_fps,dists)
    end

    #=
    #plot trajectory projection
    pl_proj = plot(traj[:, 1], traj[:, 3], ylims=ylims[1],
        xlims=xlimits_proj[i], ylabel=L"x_2", xlabel=L"x_1",
        xformatter=:auto, labels=labels, lw=2, la=1.0,label="",
        aspect_ratio=1,lc=:gray40,title=L"w=%$w_rounded",size=(600,600))
    =#

    title = titles[i]
    #plot timeseries 
    labels = [L"x_1" L"x_2" L"x_3"]
    x_vars = Matrix(traj[:][:,1:2:2N]) 
    pl_ts = plot(t, x_vars, ylims=ylims[2],lw=lw,
        xformatter=:auto,ylabel=L"x(t)",labels=labels,
        title=title, left_margin=5Plots.mm,
        tickfontsize=10)

    #plot the distance to fixed points along the trajectory
    fp_labels = [L"u_1^{*(1)}" L"u_1^{*(2)}" L"u_1^{*(3)}"]
    pl_d = plot(t, stack(dists_fps), ylims=ylims[3],
        xlabel=L"t",ylabel=L"d \left( \mathbf{u}(t),\mathbf{u}^* \right)",
        labels=fp_labels, lw=lw, tickfontsize=10,
        title=title)

    #plot trajectory projection
    pl_3d = plot(traj[:, 1], traj[:, 3], traj[:, 5],  
        ylabel=L"x_2", xlabel=L"x_1", zlabel=L"x_3",
        xformatter=:auto, labels=labels, lw=4, la=1.0, legend=false,
        aspect_ratio=1, lc=:gray40, title=title, proj_type=:ortho,
        cam=(160,30),left_margin=-25Plots.mm,
        tickfontsize=1)


    push!(plots, pl_ts, pl_d, pl_3d)

end
hline!(plots[4], [x_sw_val], lc=:gray10, lw=2 ,ls=:dash,label="")
pl = plot(plots..., layout=(2, 3), size=(1200, 600), guidefontsize=18, legendfontsize=15, dpi=300)

figs_dir = "../../paper_figs/"
savefig(pl, figs_dir * "distance_to_fixedpoints_N3_g5.png")

#plot!(pl_ts,pl_d,layout=(2,1),size=(400,600),
#    guidefontsize=18,legendfontsize=15,tickfontsize=10)

