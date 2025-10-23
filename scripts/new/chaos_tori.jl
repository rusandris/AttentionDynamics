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

ds = create_ds(u01,p;diffeq=diffeq)

Ttr = 10000
T_traj = 100
T_λ = 2000
n_psos = 5000

w_vals = [0.42,0.45,0.47,0.5]
plots = []
plots_psos = []

for w in w_vals
    @show w
    set_parameter!(ds,2,w)

    #trajectories
    traj1,t = trajectory(ds,T_traj,u01;Ttr=Ttr,Δt=1e-2)


    #psos 
    psos_points = []
    pmap = PoincareMap(ds, (1, 1.0))

    #=
    #transient
    while pmap.tcross < Ttr
        step!(pmap)
    end
    =#

    for i in 1:n_psos
        step!(pmap)
        point = current_state(pmap)[[3, 5]]
        push!(psos_points, point)
    end

    psos_pointsx = [psos_points[i][1] for i in 1:length(psos_points)]
    psos_pointsy = [psos_points[i][2] for i in 1:length(psos_points)]
    pl_psos = plot(psos_pointsx, psos_pointsy, st=:scatter,
        ms=0.5, mc=:gray10, ma=0.4, legend=false, xlabel=L"x_2")
    push!(plots_psos, pl_psos)

    #second traj (attractor)
    traj2,t = trajectory(ds,T_traj,u02;Ttr=Ttr,Δt=1e-2)

    pl = plot(traj1[:,1],traj1[:,3],traj1[:,5],label=L"1");
    plot!(traj2[:,1],traj2[:,3],traj2[:,5],label=L"2");

    plot!(pl,xlabel=L"x_1",ylabel=L"x_2",zlabel=L"x_3")
    
    #lyap exp
    λ1 = lyapunov(ds,T_λ;u0=u01,Ttr=Ttr)
    λ2 = lyapunov(ds,T_λ;u0=u02,Ttr=Ttr)
    λ1_round = round(λ1,digits=4)
    λ2_round = round(λ2,digits=4)

    #add title to plot
    plot!(camera=(130,40),legend=:bottomright,
        title=L"N=%$(N),g=%$(g),w=%$(w)"*"\n" * L"\lambda_1 = %$(λ1_round),\lambda_2 = %$(λ2_round)")

    push!(plots,pl)

end

#phase portrait axis adjustments

#psos plot adjustments
plot!(plots_psos[1],ylabel=L"x_3")
plot!(plots_psos[4],ylims=(0,3),xlims=(0,3))

pl = plot(plots...,plots_psos...,layout=(2,length(w_vals)),size=(1200,600),bottom_margin=5Plots.mm,left_margin=5Plots.mm,dpi=300)
fig_dir = "../../paper_figs/"
savefig(pl,fig_dir*"phase_portraits_3d_N=$(N)_g=$(g)_Ttraj_$(T_traj)_T_lambda_$(T_λ)_TTr_$(Ttr)"*".png")


#traj distance scaling with time
#test for partially predictable chaos?

#=
predictability(ds; T_max=2000,n_samples=30,Ttr = 500,T_sample=1e4)

#w=w_vals[1]
#(:SC, -0.0023257727854968857, -0.08707304141994654)

#w=w_vals[2]
#(:REG, 0.9537348139942716, 0.9999999999981031)

#w=w_vals[3]
#(:PPC, 0.15868467475768436, 0.9308109639269125)

#w=w_vals[4]
#(:REG, 1.0849091176297923, 0.9999999991445528)
=#

plot_d = plot()
for w in w_vals
    @show w
    set_parameter!(ds,2,w)

    #traj distances 
    traj_dists1,t = traj_dist(ds,u01,5000;Ttr=Ttr,d0=1e-12/sqrt(2N))
    #traj_dists2 = traj_dist(ds,u02,T_traj)
    plot!(plot_d,t,traj_dists1,yaxis=:log,label="w=$w")
    
end

plot!(plot_d,xlabel=L"t",ylabel=L"d(\mathbf{u})",guidefontsize=15,legend=:bottomright,dpi=300)
savefig(plot_d,fig_dir*"trajectory_distance_scaling_3d_N=$(N)_g=$(g)_Ttraj_$(T_traj)_T_lambda_$(T_λ)_TTr_$(Ttr)"*".png")