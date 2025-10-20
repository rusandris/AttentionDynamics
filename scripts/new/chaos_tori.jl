using Revise
using BifurcationKit
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using Plots,LaTeXStrings
using DynamicalSystems
cd(@__DIR__)

#-----------IDEA---------------
#plot torus-chaos phases for several dimensions 
#lyap exponent in time? trajectory distance in time

function traj_dist(ds,u0,T;Ttr=1000,d0 = 1e-9 / sqrt(dimension(ds)),Δt=1e-2)

    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ rand(dimension(ds)) .* d0])

    #transient
    while current_time(pds) < Ttr
        step!(pds,Δt)
    end

    #rescale 
    u = current_states(pds)[1]
    pds = ParallelDynamicalSystem(ds,[u,u .+ rand(dimension(ds)) .* d0])

    #distance
    t = Float64[]
    dists = Float64[]
    while current_time(pds) < T
        step!(pds,Δt)
        u0,up = current_states(pds)
        d = norm(u0-up)
        push!(dists,d)
        push!(t,current_time(pds))
    end

    return dists,t
end


N = 6
g = 12.0

wc1 = (N-1)/g
wc2 = (N-1)/N
w = 0.5
p = [g,w,N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
u01 = rand(2N)
u02 = rand(2N)

ds = create_ds(u0,p;diffeq=diffeq)

Ttr = 10000
T_traj = 100
T_λ = 2000

w_vals = [0.42,0.45,0.47,0.5]
plots = []
plot_d = plot()
for w in w_vals
    @show w
    set_parameter!(ds,2,w)

    #psos 
    #=
        psos_points = []
        pmap = PoincareMap(ds, (1,1.0))
        psos_pointsx = [psos_points[i][1] for i in 1:length(psos_points)]
        psos_pointsy = [psos_points[i][2] for i in 1:length(psos_points)]
        plot(psos_pointsx[1:150],psos_pointsy[1:150],st=:scatter,ms=0.5,mc=:gray10,ma=0.4)
    =#

    #=
    #trajectories
    traj1,t = trajectory(ds,T_traj,u01;Ttr=Ttr,Δt=1e-2)
    traj2,t = trajectory(ds,T_traj,u02;Ttr=Ttr,Δt=1e-2)

    pl = plot(traj1[:,1],traj1[:,3],traj1[:,5],label=L"1");
    plot!(traj2[:,1],traj2[:,3],traj2[:,5],label=L"2");

    
    #lyap exp
    λ1 = lyapunov(ds,T_λ;u0=u01,Ttr=Ttr)
    λ2 = lyapunov(ds,T_λ;u0=u02,Ttr=Ttr)
    λ1_round = round(λ1,digits=4)
    λ2_round = round(λ2,digits=4)

    #add title to plot
    plot!(camera=(130,40),legend=:bottomright,
        title=L"N=%$(N),g=%$(g),w=%$(w)"*"\n" * L"\lambda_1 = %$(λ1_round),\lambda_2 = %$(λ2_round)")

    push!(plots,pl)
    =#
    #traj distances 
    traj_dists1,t = traj_dist(ds,u01,5000;Ttr=Ttr,d0=1e-12/sqrt(2N))
    #traj_dists2 = traj_dist(ds,u02,T_traj)
    plot!(plot_d,t,traj_dists1,yaxis=:log,label="w=$w")

end
pl = plot(plots...,layout=(1,length(w_vals)),size=(1200,350),dpi=300)
fig_dir = "../../paper_figs/"
savefig(pl,fig_dir*"phase_portraits_3d_N=$(N)_g=$(g)_Ttraj_$(T_traj)_T_lambda_$(T_λ)_TTr_$(Ttr)"*".png")

plot!(plot_d,xlabel=L"t",ylabel=L"d(\mathbf{u})",guidefontsize=15,legend=:bottomright,dpi=300)
savefig(plot_d,fig_dir*"trajectory_distance_scaling_3d_N=$(N)_g=$(g)_Ttraj_$(T_traj)_T_lambda_$(T_λ)_TTr_$(Ttr)"*".png")