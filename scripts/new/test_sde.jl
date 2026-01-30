using Revise
using AttentionDynamics
using StochasticDiffEq
using OrdinaryDiffEq
using Plots,LaTeXStrings
using DynamicalSystems
using DiffEqNoiseProcess
cd(@__DIR__)

#diagonal noise function
function noise_function!(du, u, p, t) 
    σ = p[4]
    N = Int(p[3])
    for i in 1:2:(2*N)
        du[i] = σ #*u[i]
    end
end

function noise_function_orig!(du, u, p, t)
    σ = p[6]
    N = Int(p[3])
    for i in 1:2:(2*N)
        du[i] = σ*u[i]
    end
end

function noise(u, p, t)
    return abs(rand())
end

noise_process = NoiseFunction(0.0, noise)

N = 10
g = 18.0
w = 0.7
σ = 1e-8 #5*1e-1 #5*1e-3
u0 = rand(2N) .* 0.1
p = [g,w,N,σ]
T = 1800
tspan = (0.0,T)
xlims=(0.0,180.0)

σs = [1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
σs = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
diffeq = (alg=Tsit5(), dtmax=1e-2, adaptive=true, maxiters=1e15)
plots = []
for σ in σs
    pl = plot(ylabel=L"x_i", xlabel="t", guidefontsize=12,
        legendfontsize=5)

    #noiseless rescaled 
    ds = create_ds(u0, [g, w, N]; diffeq=diffeq) #no noise
    traj, t = trajectory(ds, T; Ttr=0, Δt=1e-2)
    plot!(t, traj[:, 1], lc=:gray10, lw=1, xlims=xlims,
        label=L"\sigma=0")


    #rescaled equations as SDE
    prob_sde = SDEProblem(att_din!, noise_function!, u0, tspan,p;noise=noise_process)
    sol = solve(prob_sde;maxiters=1e7,dtmax=1e-2)
    @show sol.retcode
    #plot(sol,legend=false)
    #x1 = [sol.u[i][1] for i in 1:length(sol.u)] 
    #plot!(sol.t, x1, lw=1, xlims=xlims,lc=:blue,
    #    label="SDE, " * L"\sigma=%$σ")


    #noisey rescaled 
    diffeq = (alg=Tsit5(), dtmax=1e-2, adaptive=true, maxiters=1e15)
    ds_noise = create_ds_noisy(u0, [g, w, N, σ]; diffeq=diffeq) #no noise
    traj_noise, t = trajectory(ds_noise, T; Ttr=0, Δt=1e-2)
    plot!(t, traj_noise[:, 1], lc=:green, lw=1, xlims=xlims,
        label="ODE, " * L"\sigma=%$σ")

    #noisey rescaled proportional
    diffeq = (alg=Tsit5(), dtmax=1e-2, adaptive=true, maxiters=1e15)
    ds_noise_prop = create_ds_prop_noisy(u0, [g, w, N, σ]; diffeq=diffeq) #no noise
    traj_noise, t = trajectory(ds_noise_prop, T; Ttr=0, Δt=1e-2)
    plot!(t, traj_noise[:, 1], lc=:orange, lw=1, xlims=xlims,
        label="ODEm, " * L"\sigma=%$σ")
    push!(plots,pl)
end
pl = plot(plots..., layout=(6, 1), size=(600, 1200), 
    plot_title=L"N=%$(N),g=%$(g),w=%$(w)",dpi=300,
    legend=:topright,legend_columns=-1,xlims=(1700,1800))
savefig(pl,"../../paper_figs/timeseries_transient_chaotic_noise_tests_end_long.png")    
#=
#original equations as SDE

#original params
N = 300
r_const = 9.0
rc = r_const
rp = r_const #9.0;
K = 1.0;
c = 2.4;
α = 0.005;
w, g = rescaled_params(rp, rc, K, c, α, N)
p_orig = [rp, rc, K, c, α]

wc1 = (N - 1) / (rp / α)
wc2 = (N - 1) / N

σ = 1e-7
p_orig_sde = [p_orig...,σ]
tspan = (0.0,100.0)
u0 = rand(2N) .* 0.1
prob_orig_sde = SDEProblem(att_din_original!, noise_function_orig!, u0, 
    tspan, p_orig_sde; noise=noise_process)

sol = solve(prob_orig_sde; maxiters=1e15)
x1 = [sol.u[i][1] for i in 1:length(sol.u)]
plot(sol.t, x1, lw=2)

#original 

ds_orig = create_original_ds(u0, p_orig; diffeq=diffeq)
traj, t = trajectory(ds_orig, 100; Ttr=0, Δt=1e-2)
plot(t, traj[:, 1], lc=:gray10, lw=1)
=#