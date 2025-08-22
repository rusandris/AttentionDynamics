using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
cd(@__DIR__)

figs_dir = "../../paper_figs/"

#-----------------------------N=2-------------------------------

N = 2 #N = 5
g = 5.0 #g = 20.0
w = 0.3 #w = 0.25
σ = 1e-1
p_noisy = [g, w, N, σ]
p = [g, w, N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)
ds_noisy = create_ds_noisy(u0, p_noisy; diffeq=diffeq)

σs = [1e-3, 1e-2, 1e-1, 0.0] #noise levels
colors = [:gray70, :gray40, :gray10, :red]
pl_N2 = plot(; guidefontsize=20, legendfontsize=15, xlabel=L"x_1", ylabel=L"y_1",
    title=L"N = %$N" * ", " * L"g = %$g" * L"w = %$w" * ", with added noise",
    legendtitle="uniform noise, " * L"\sigma" * " amplitude", dpi=300)

crs_means = []
crs_extrema = []
for (i, σ) in enumerate(σs)
    @show i
    set_parameter!(ds_noisy, 4, σ)
    traj_noisy, t = trajectory(ds_noisy, 150, u0; Ttr=2000, Δt=1e-3)
    plot!(traj_noisy[:, 1], traj_noisy[:, 2], lc=colors[i], label=L"\sigma = %$σ")

    #calc contraction data 
    crs = [contraction_rate(traj_noisy[i], p) for i in 1:length(traj_noisy)]
    push!(crs_means,mean(crs))
    push!(crs_extrema, extrema(crs))

end

#-----------------------------N=10-------------------------------

N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.4 #w = 0.25
σ = 1e-1
p_noisy = [g, w, N, σ]
p = [g, w, N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)
ds_noisy = create_ds_noisy(u0, p_noisy; diffeq=diffeq)

traj_noisy, t = trajectory(ds_noisy, 150, u0; Ttr=2000, Δt=1e-3)
traj, t = trajectory(ds, 50, u0; Ttr=2000, Δt=1e-3)


#plot!(traj_noisy[:, 1], traj_noisy[:, 2], lc=:gray10, label=L"\sigma = 0.1")
#plot!(traj[:, 1], traj[:, 2], lc=:red, lw=1, label=L"\sigma = 0")

#lyapunov(ds_noisy,200;Ttr=200)

σs = [1e-3,1e-2,1e-1,0.0] #noise levels
colors = [:gray70, :gray40, :gray10,:red]
pl_N10 = plot(; guidefontsize=20, legendfontsize=15, xlabel=L"x_1", ylabel=L"y_1",
    title=L"N = %$N" * ", " * L"g = %$g" * L"w = %$w" * ", with added noise", dpi=300)


crs_means = []
crs_extrema = []
for (i,σ) in enumerate(σs)
    set_parameter!(ds_noisy,4,σ)
    traj_noisy, t = trajectory(ds_noisy, 150, u0; Ttr=2000, Δt=1e-3)
    plot!(traj_noisy[:, 1], traj_noisy[:, 2], lc=colors[i], label=L"\sigma = %$σ")    


    #calc contraction data 
    crs = [contraction_rate(traj_noisy[i], p) for i in 1:length(traj_noisy)]
    push!(crs_means, mean(crs))
    push!(crs_extrema, extrema(crs))
end
pl = plot(pl_N2, pl_N10, size=(800, 300), left_margin=5Plots.mm, 
    bottom_margin=5Plots.mm, legendfontsize=8, legendtitlefontsize=6, 
    legend_columns=1, legend=:bottomright,titlefontsize=10)
savefig(pl,figs_dir*"noisy_equations_test.png")