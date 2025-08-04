using Plots,LaTeXStrings
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
#using Attractors
cd(@__DIR__)
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

#traj
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.3
p = [g,w,N]

T_plot_chaos = 400
Ttr_chaos = 4000
Δt = 0.005
Tλ = 200

pl = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
        guidefontsize=20,legend=:bottomright,dpi=300)

pl_cr = plot(xlabel=L"t",ylabel=L"Trace(J)",
        guidefontsize=20,legend=:bottomright,dpi=300)

rand_amplitudes = [1,50]

for j in 1:length(rand_amplitudes) 
    A = rand_amplitudes[j]

    for i in 1:5

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

        plot!(pl,x1 .- x2, y1 .- y2,
                lw=1, la=0.7, c=cgrad(:default), clims=(0.0, 30.0), line_z=λ,
                ylims=(-4,4),xlims=(-25,25),label=L"A = %$A")


        #contraction (trace of Jacobian along trajectory)
        trJ = [contraction_rate(traj_chaos[i],p) for i in 1:length(traj_chaos)]
        @show extrema(trJ)
        @show mean(trJ)

        plot!(pl_cr,t, trJ,
            lw=1, la=1.0, c=cgrad(:default), clims=(0.0, 30.0), line_z=λ,
            ylims=(0,50),xlims=(100,110),label=L"A = %$A")
    end
end
pl
pl_cr
savefig(pl, figs_dir*"chaotic_regime_diff_proj_lyap_colored_test_Tlamb$(Tλ)_Ttr$(Ttr_chaos)_N$N"*"_g$g"*"_w$w"*".png")
#IDEA: compute unstable limit cycles for this parameter?
#using Bifkit.jl? 