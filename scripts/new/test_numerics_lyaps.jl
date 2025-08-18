using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
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
Ttr = 400
nr_u0s = 10
d0 = 1e-9

#def initial conditions
u0s = [rand(2N) .* 50 for i in 1:20]
lyaps = []

#----------------------------define ds------------------------------- 
begin 
    diffeq0 = (alg = Tsit5(),dt = 1e-3,adaptive=true,maxiters=1e15) 
    diffeq1 = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15) 
    diffeq2 = (alg = Rodas5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15) 
    diffeq2_hp = (alg = Rodas5(),dt = 1e-4,dtmax=1e-4,adaptive=true,maxiters=1e15) 
end
diffeqs = [diffeq0,diffeq1,diffeq2,diffeq2_hp]

systems = []
for diffeq in diffeqs 
    ds = create_ds(u0s[1],p;diffeq=diffeq)
    push!(systems,ds)
end

#-------------------------tests-----------------------
T_short = 5
T_limcycle = 20
ylims = (-3,3)

titles = ["Tsit5, adaptive","Tsit5, adaptive, dtmax=1e-3",
    "Rodas5, adaptive, dtmax=1e-3","Rodas5, adaptive, dtmax=1e-4"]
plots = []

for (i,ds) in enumerate(systems[end:end])
    @show i
    traj,t = trajectory(ds,T_limcycle,u0s[1];Ttr=Ttr,Δt = 0.001)
    pl = plot_projection(traj;xlims=(-30,30),ylims=ylims,lw=4,lc=:gray50,
        label="partial limit cycle",title= titles[i], #L"N=%$(N)"* L",g=%$(g)" *
        xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",guidefontsize=20,legendfontsize=6) 

    println("1")
    #see which variables are nonzero
    active_vars_search = findall(u -> u > 1.0, maximum(Matrix(traj);dims=1))
    active_vars = [last(ind.I) for ind in active_vars_search]

    #perturb trajectory on the limit cycle in the whole space
    u_p = inittest_custom(traj[end],d0;dims=[1:20;])
    traj_p,t = trajectory(ds,T_short,u_p;Ttr=0,Δt = 0.001)
    plot_projection!(pl,traj_p;xlims=(-30,30),ylims=ylims,label="perturbed uniformly")
    println("1")

    #perturb trajectory on the limit cycle only in the invariant subspace
    u_p_inv = inittest_custom(traj[end],d0;dims=active_vars)
    traj_p_inv,t = trajectory(ds,T_short,u_p_inv;Ttr=0,Δt = 0.001)
    plot_projection!(pl,traj_p_inv;xlims=(-30,30),ylims=ylims,label="perturbed in invariant subspace")
    println("1")

    push!(plots,pl)
end

plot(plots...,layout=(length(plots),1),size=(500,1200))

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