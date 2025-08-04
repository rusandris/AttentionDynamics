using Plots,LaTeXStrings
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
#using Attractors
cd(@__DIR__)
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

#traj
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.25
p = [g,w,N]

T_plot_chaos = 400
Ttr_chaos = 14000
Δt = 0.005
Tλ = 200

pl = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
        guidefontsize=20,legend=:bottomright,dpi=300)

for i in 1:10

        u0 = rand(2N) #.* 50
        ds = create_ds(u0,p;diffeq=diffeq)

        #plot transient
        traj_chaos,t_tr = trajectory(ds,T_plot_chaos,u0;Δt = Δt)

        #lyapunov exponent differs?
        λ = lyapunov(ds,Tλ;u0=current_state(ds))
        @show λ

        x1 = traj_chaos[:,1]
        x2 = traj_chaos[:,3]
        y1 = traj_chaos[:,2]
        y2 = traj_chaos[:,4]

        plot!(x1 .- x2, y1 .- y2,
                lw=1, la=0.7, c=cgrad(:default), clims=(0.0, 10.0), line_z=λ,
                label="",ylims=(-4,4),xlims=(-25,25))
end
pl
savefig(pl, figs_dir*"chaotic_regime_diff_proj_scattered_test_N$N"*"_g$g"*"_w$w"*".png")
#IDEA: compute unstable limit cycles for this parameter?
#using Bifkit.jl? 