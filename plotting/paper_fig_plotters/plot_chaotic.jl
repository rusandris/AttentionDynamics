using Plots,LaTeXStrings
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
#using Attractors
cd(@__DIR__)
figs_dir = "../../paper_figs/"
mkpath(figs_dir)


function plot_projection(traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf))
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    pl = plot(x1 .- x2, y1 .- y2,
        lw=1, la=1.0, ylims=ylims, xlims=xlims)

end

function plot_projection!(pl,traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf))
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    plot!(pl,x1 .- x2, y1 .- y2,
        lw=1, la=1.0, ylims=ylims, xlims=xlims)

end

#traj
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.3
p = [g,w,N]

T_plot_chaos = 400
Ttr_chaos = 7000
Δt = 0.005
Tλ = 200

pl = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
        guidefontsize=20,legend=:bottomright,dpi=300)

pl_cr = plot(xlabel=L"t",ylabel=L"Trace(J)",
        guidefontsize=20,legend=:bottomright,dpi=300)

rand_amplitudes = [1,50]
nr_ic = 10
plots_all = []
for j in 1:length(rand_amplitudes) 
    plots = []
    A = rand_amplitudes[j]
    for i in 1:nr_ic

        u0 = rand(2N) .* A
        ds = create_ds(u0,p;diffeq=diffeq)

        #plot transient
        traj_chaos,t = trajectory(ds,T_plot_chaos,u0;Ttr = Ttr_chaos,Δt = Δt)

        #lyapunov exponent differs?
        λ = lyapunov(ds,Tλ;u0=current_state(ds))
        @show λ

        x1 = traj_chaos[:,1]
        x2 = traj_chaos[:,3]
        y1 = traj_chaos[:,2]
        y2 = traj_chaos[:,4]

        pl = plot(x1 .- x2, y1 .- y2,
            lw=1, la=1.0, ylims=(-4, 4), xlims=(-25, 25),
                label=L"A = %$A")
        
        #c=cgrad(:default), clims=(0.0, 30.0), line_z=λ
        push!(plots,pl)
        #=
        #contraction (trace of Jacobian along trajectory)
        trJ = [contraction_rate(traj_chaos[i],p) for i in 1:length(traj_chaos)]
        @show extrema(trJ)
        @show mean(trJ)

        plot!(pl_cr,t, trJ,
                lw=1, la=1.0, c=cgrad(:default), clims=(0.0, 30.0), line_z=λ,
                ylims=(0,50),xlims=(100,110),label=L"A = %$A")
        =#
    end
   pl =  plot(plots..., layout=(2,5),size=(1600,400),dpi=300)
   push!(plots_all,pl)
end
savefig(plots_all[2],figs_dir*"chaotic_scattered_test.png")
pl_cr
#savefig(pl, figs_dir*"chaotic_regime_diff_proj_lyap_colored_test_Tlamb$(Tλ)_Ttr$(Ttr_chaos)_N$N"*"_g$g"*"_w$w"*".png")
#IDEA: compute unstable limit cycles for this parameter?
#using Bifkit.jl? 


#test partial limit cycle trajectories
#perturbed inside invariant plane vs. outside
#=
function inittest_custom(u1,d0;dims=[3,5,9])
           u_p = zeros(length(u1))
           u_p[dims] .+= rand(length(dims)) .* d0 
           return u1 .+ u_p
       end

u0 = rand(2N) .* 50
traj,t = trajectory(ds,500,u0;Ttr=1000,Δt = 0.001)
M = Matrix(traj)
minmaxs = [extrema(M[:,i]) for i in 1:2:20]
d0 = 1e-9
u0_limc = traj[end]
u0_p_inv = inittest_custom(traj[end],d0)
u0_p = inittest_custom(traj[end],d0;dims=[1:2N;])
   
traj_limc,t = trajectory(ds,5,u0_;Ttr=0,Δt = 0.001)
traj_p,t = trajectory(ds,5,u0_p;Ttr=0,Δt = 0.001)
traj_p_inv,t = trajectory(ds,5,u0_p_inv;Ttr=0,Δt = 0.001)

pl = plot_projection(traj_limc;xlims=(-30,30),ylims=(-4,4))
plot_projection!(pl,traj_p;xlims=(-30,30),ylims=(-4,4))
plot_projection!(pl,traj_p_inv;xlims=(-30,30),ylims=(-4,4))

#check Jacobian along limit cycle
J = build_jacobian(u0_limc,p)
vals,vecs = eigen(J)

=#