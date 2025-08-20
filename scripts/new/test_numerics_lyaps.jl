using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DelimitedFiles
using DynamicalSystems
cd(@__DIR__)
#inspectdr()
gr()

function plot_projection(traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf),kwargs...)
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    pl = plot(x1 .- x2, y1 .- y2;lw=1, la=1.0, 
        ylims=ylims, xlims=xlims,kwargs...)

end

function plot_projection!(pl,traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf),kwargs...)
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    plot!(pl,x1 .- x2, y1 .- y2; lw=1, la=1.0, 
        ylims=ylims, xlims=xlims,kwargs...)

end

function inittest_custom(u1,d0=1e-9;dims=[3,4,5,6,9,10])
    D = length(u1)
    u_p = zeros(D)
    u_p[dims] .+= rand(length(dims)) .* d0/sqrt(D) 
    return u1 .+ u_p
end

#check if partial limit cycles really exist


#---------------------system params----------------------------
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.3 #w = 0.25

p = [g, w, N]

#search a u0 before transition to limit cycle
T = 50
Ttr = 50
nr_u0s = 10
d0 = 1e-9
maxiters = 1e8

#def initial conditions
#u0s = [rand(2N) .* 50 for i in 1:20]
u0 = vec(readdlm("../../data/initial_condition_N10.txt"))
u0s = [u0]
lyaps = []

#----------------------------define ds------------------------------- 
begin 
    #Tsit5
    diffeq0 = (alg = Tsit5(),adaptive=true,maxiters=maxiters) 
    diffeq1 = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=maxiters) 
    diffeq2 = (alg = Tsit5(),dt = 1e-3,adaptive=false,maxiters=maxiters)
    diffeq3 = (alg = Tsit5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=maxiters) 
    diffeq4 = (alg = Tsit5(),dt = 1e-4,adaptive=false,maxiters=maxiters) 

    #Rodas5
    diffeq5 = (alg = Rodas5(),adaptive=true,maxiters=maxiters)
    diffeq6 = (alg = Rodas5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=maxiters)
    diffeq7 = (alg = Rodas5(),dt = 1e-3,adaptive=false,maxiters=maxiters)
    diffeq8 = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=maxiters)
    diffeq9 = (alg = Rodas5(),dt = 1e-4,adaptive=false,maxiters=maxiters)
end
diffeqs = [diffeq0,diffeq1,diffeq2,diffeq3,diffeq4,
    diffeq5,diffeq6,diffeq7,diffeq8,diffeq9]

systems = []
for diffeq in diffeqs 
    ds = create_ds(u0s[1],p;diffeq=diffeq)
    push!(systems,ds)
end


#-------------------------tests-----------------------
T_short = 5
T_limcycle = 20
ylims = (-3,3)
lws = [1,5,3,2,1]

titles = ["adaptive, dtmax = Inf","adaptive, dtmax=1e-3",
    "fixed, dt=1e-3","adaptive, dtmax=1e-4","fixed, dt=1e-4"]

pl_Tsit5 = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
    guidefontsize=20,legendfontsize=6,title="N=$(N), g=$(g), Tsit5")
pl_Rodas5 = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
    guidefontsize=20,legendfontsize=6,title="N=$(N), g=$(g), Rodas5")

colors = distinguishable_colors(5)

#plot test trajectories first 
for (i,ds) in enumerate(systems)
    @show i
    i == 6 && continue

    traj,t = trajectory(ds,T_limcycle,u0s[1];Ttr=Ttr,Δt = 0.001)

    if i < 6
        plot_projection!(pl_Tsit5,traj;xlims=(-30,30),ylims=ylims,
        label = titles[i],lw=lws[i],lc=colors[i])
    else
        plot_projection!(pl_Rodas5,traj;xlims=(-30,30),ylims=ylims,
        label = titles[i-5], lw=lws[i-5],lc=colors[i-5])
    end
end
plot!(pl_Tsit5,legend=:bottomright)
plot!(pl_Rodas5,legend=:bottomright)
pl = plot(pl_Tsit5,pl_Rodas5,size=(600,400),dpi=300)
savefig(pl,figs_dir*"solver_comparison_Tsit5_vs_Rodas5_test.png")

perturb_plots = []
clims = (0.0,2.0)

for (i,ds) in enumerate(systems)

    @show i

    
    i == 6 && continue
    if i < 6
        plot_title = "Tsit5, " * titles[i]
    else 
        plot_title = "Rodas5, " * titles[i-5]
    end

    traj,t = trajectory(ds,T_limcycle,u0s[1];Ttr=Ttr,Δt = 0.001)
    pl = plot_projection(traj;xlims=(-30,30),ylims=ylims,lw=3,lc=:gray70,
        label="partial limit cycle",title = plot_title, #L"N=%$(N)"* L",g=%$(g)" *
        xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",guidefontsize=20,legendfontsize=6) 

    #see which variables are nonzero
    active_vars_search = findall(u -> u > 1.0, maximum(Matrix(traj);dims=1))
    active_vars = [last(ind.I) for ind in active_vars_search]

    #perturb trajectory on the limit cycle in the whole space
    u_p = inittest_custom(traj[end],d0;dims=[1:20;])
    λ_p =  lyapunov(ds,100;u0=traj[end],inittest = (u1,d0) -> u_p,Ttr=100)
    @show λ_p
    traj_p,t = trajectory(ds,T_short,u_p;Ttr=0,Δt = 0.001)
    plot_projection!(pl,traj_p;xlims=(-30,30),ylims=ylims,label="perturbed uniformly",
        c=cgrad(:Set1_3), clims=clims, line_z=λ_p)


    #perturb trajectory on the limit cycle only in the invariant subspace
    u_p_inv = inittest_custom(traj[end],d0;dims=active_vars)
    traj_p_inv,t = trajectory(ds,T_short,u_p_inv;Ttr=0,Δt = 0.001)
    λ_p_inv =  lyapunov(ds,100;u0=traj[end],inittest = (u1,d0) -> u_p_inv,Ttr=100)
    @show λ_p_inv
    plot_projection!(pl,traj_p_inv;xlims=(-30,30),ylims=ylims,label="perturbed in invariant subspace",
        c=cgrad(:Set1_3), clims=clims, line_z=λ_p_inv)


    push!(perturb_plots,pl)
end

figs_dir = "../../numerics_figs/"
mkpath(figs_dir)
#plot(perturb_plots...,layout=(5,5),size=(200,300))
for (i,diffeqs) in enumerate(perturb_plots)
    i == 6 && continue
    pl = perturb_plots[i]
    plot!(pl,dpi=300)
    savefig(pl,figs_dir*"perturb_trajectory_tests_pl$i"*".png")

end


#calc lyapunov exponents using different perturbation vectors
#lyapunov(ds,100;u0=u0_limc,inittest = (u1,d0) -> inittest_custom(u1,d0;dims=[1,3,9,17]),Ttr=100)
#lyapunov(ds,100;u0=u0_limc,inittest = (u1,d0) -> inittest_custom(u1,d0;dims=[1:2N;]),Ttr=100)


#check Jacobian along limit cycle
#J = build_jacobian(u0_limc,p)
#vals,vecs = eigen(J)


#=
for u0 in u0s[1:2]

    λ1 = lyapunov(ds1, 100; u0=u0, Ttr=Ttr)
    @show λ1
    λ2 = lyapunov(ds2, 100; u0=u0, Ttr=Ttr)
    @show λ2
    λ2_hp = lyapunov(ds2_hp, 100; u0=u0, Ttr=Ttr)
    @show λ2_hp
    push!(lyaps,λ1,λ2,λ2_hp)

end
=#