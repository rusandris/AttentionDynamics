using Revise
using Plots, LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
cd(@__DIR__)

figs_dir = "../../paper_figs/"

T = 100 
Ttr = 5000

#-----------------------------N=2-------------------------------

#system setup
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
σs = [0.0, 1e-3, 5 * 1e-3, 1e-2, 5 * 1e-2, 1e-1] #[1e-3, 1e-2, 1e-1, 0.0] #[0.0,1e-3,5*1e-2]  #noise levels
trajN2, t = trajectory(ds, 0, u0; Ttr=Ttr, Δt=1e-3) #start other trajs from here
u0_traj = trajN2[end] #on the attractor

#-------------------------plot------------------------
#colors = distinguishable_colors(length(σs);lchoices=range(30,50;length=15)) #[:gray80, :gray60, :gray40, :gray20, :gray10, :red]
colors = range(colorant"gray10", stop=colorant"orange", length=18)[1:6:end] 
pl_N2 = plot(; guidefontsize=20, legendfontsize=15, 
    xlabel=L"x_1", ylabel=L"h_1",dpi=300,legend=:topright,
    xlims = (0.0,4.0),ylims = (0.0,4.0));
#legendtitle="uniform noise, " * L"\sigma" * " amplitude"

crs_means = []
crs_extrema = []
@info "N=2 system (no chaos)"
for (i, σ) in enumerate(σs)
    @show σ
    set_parameter!(ds_noisy, 4, σ)
    traj_noisy, t = trajectory(ds_noisy, T, u0_traj; Ttr=Ttr, Δt=1e-3)
    plot!(traj_noisy[:, 1], traj_noisy[:, 2], lc=colors[i], 
        label=L"\sigma = %$σ",lw=2)

    #calc contraction data 
    #crs = [contraction_rate(traj_noisy[i], p) for i in 1:length(traj_noisy)]
    #push!(crs_means,mean(crs))
    #push!(crs_extrema, extrema(crs))
end
pl_N2
#-----------------------------N=10-------------------------------

#system setup
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.4 #w = 0.25
σ = 1e-1
p_noisy = [g, w, N, σ]
p = [g, w, N]

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)
ds = create_ds(u0, p; diffeq=diffeq)

#starting init cond
trajN2, t = trajectory(ds, 0, u0; Ttr=Ttr, Δt=1e-3) #start other trajs from here
u0_traj = trajN2[end] #on the attractor

ds_noisy = create_ds_noisy(u0, p_noisy; diffeq=diffeq)

#-------------------------plot------------------------
pl_N10 = plot(; guidefontsize=20, legendfontsize=15, xlabel=L"x_1", 
    ylabel=L"h_1", dpi=300,legend=:topright,yticks=[0,1,2,3,4],
    ylims=[0.0,4.0])


crs_means = []
crs_extrema = []
@info "N=10 system (chaos)"
for (i,σ) in enumerate(σs)
    @show σ
    set_parameter!(ds_noisy,4,σ)
    traj_noisy, t = trajectory(ds_noisy, T, u0_traj; Ttr=Ttr, Δt=1e-3)
    plot!(traj_noisy[:, 1], traj_noisy[:, 2], lc=colors[i],
        label=L"\sigma = %$σ",lw=2)    

    #calc contraction data 
    #crs = [contraction_rate(traj_noisy[i], p) for i in 1:length(traj_noisy)]
    #push!(crs_means, mean(crs))
    #push!(crs_extrema, extrema(crs))
end

#add fig label annotations
annotate!(pl_N2, -0.5, 3.8, text(L"(a)", :gray10, :center, 20))
annotate!(pl_N10, -5, 3.8, text(L"(b)", :gray10, :center, 20))

pl = plot(pl_N2, pl_N10, size=(800, 300), left_margin=5Plots.mm, 
    bottom_margin=5Plots.mm, legendfontsize=10, 
    legend_columns=3)

#savefig(pl,figs_dir*"noisy_equations_same_initcond.png")
savefig(pl, figs_dir * "noisy_equations_T$(T)_Ttr$(Ttr)_dt=1e-3_initcond_on_attr_black_orange.png")