using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using OrdinaryDiffEq
cd(@__DIR__)

#original params
N = 10
g = 45.0
rc = 1.0;K = 1.0;α = 5*1e-2
competition_strength(w,rc,K,α,N) = w/(1-w) * rc/(α*K*(N-1))
c = competition_strength(0.3, rc, K, α, N)

production_rate(g,α) = g*α
rp = production_rate(g, α)
p_orig = [rp,rc,K,c,α]

#=
#plot param conversion function
rcs = [1:1:100;];rps = [0.0:0.1:1;];cs = [1:0.1:10;];αs = [1e-2:1e-3:1e-1;]
pl1 = plot(rcs,first.(rescaled_params.(rp,rcs,K,c,α,N)),xlabel=L"r_c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl2 = plot(rps,last.(rescaled_params.(rps,rc,K,c,α,N)),xlabel=L"r_p",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl3 = plot(cs,first.(rescaled_params.(rp,rc,K,cs,α,N)),xlabel=L"c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl4 = plot(αs,first.(rescaled_params.(rp,rc,K,cs,α,N)),xlabel=L"\alpha",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl5 = plot(αs,last.(rescaled_params.(rp,rc,K,c,αs,N)),xlabel=L"\alpha",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl_param_conversion = plot(pl1,pl3,pl4,pl5,pl2,layout=(5,1),size=(500,1500),dpi=300,left_margin=8Plots.mm)
figs_dir = "../../paper_figs/"
savefig(pl_param_conversion,figs_dir*"pl_param_conversion_test.png")
=#

#rescaled params
w,g = rescaled_params(rp,rc,K,c,α,N)
f1, f2 = rescaling_factors(rp, rc, K, c, α, N)
g_normed = g/(N-1)
p = [g,w,10]

#initialize systems
u0 = rand(2N) .* 50
U0 = u0 .* f1
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
diffeq2 = (alg=Rodas5(), dt=1e-4, dtmax=1e-4, adaptive=true, maxiters=1e8)

#diffeq = (alg=Tsit5(), reltol=1e-8, abstol=1e-10)
ds = create_ds(u0,p;diffeq=diffeq)
ds2 = create_ds(u0, p; diffeq=diffeq2)
ds_orig = create_original_ds(U0, p_orig; diffeq=diffeq)

@time traj,t = trajectory(ds,50,u0;Ttr=300,Δt=0.001)
@time traj2, t = trajectory(ds2, 50, u0; Ttr=300, Δt=0.001)
@time traj_orig, t_orig = trajectory(ds_orig, 1000, U0; Ttr=6000, Δt=0.01)

#timeseries
pl = plot(t,traj[:,1],lc=:gray10,xlabel=L"t",ylabel=L"x_1",label="Tsit5, dt=1e-3",guidefontsize=20)
pl2 = plot(t, traj2[:, 1], lc=:gray10, xlabel=L"t", ylabel=L"x_1", label="Rodas5, dt=1e-4", guidefontsize=20)

pl_orig = plot(t_orig,traj_orig[:,1],lc=:gray10,xlabel=L"t",ylabel=L"L_1",label="Tsit5, dt=1e-3",guidefontsize=20)
pl_ts = plot(pl,pl2,pl_orig;layout=(3,1),size=(800,400),left_margin=5Plots.mm,
    plot_title = L"r_p = %$(rp), r_c = %$(rc), K = %$(K),c \approx %$(round(c;digits=3)),\alpha = %$(α)"  )

figs_dir = "../../paper_figs/"

savefig(pl_ts,figs_dir*"orig_vs_modelled_ts_c_$c"*".png")

#phase space
pp = plot(traj[:, 1], traj[:, 2], xlabel=L"x_1", 
    ylabel=L"y_1", label="Tsit5, dt=1e-3,Ttr = 300", lc=:red,
    guidefontsize=20, title=L"g = %$(g), w = %$(round(w;digits=1))")

plot!(pp,traj2[:, 1], traj2[:, 2], lc=:gray10, xlabel=L"x_1", 
    ylabel=L"y_1", label="Rodas5, dt=1e-4, Ttr = 300",
    guidefontsize=20,
    ylims=(0,2.5))

pp_orig = plot(traj_orig[:, 1], traj_orig[:, 2], lc=:gray10, xlabel=L"L_1", 
    ylabel=L"Y_1", label="Tsit5, dt=1e-3, Ttr = 6000", guidefontsize=20, ylims=(0, 2.5),legend=:bottomright,
    title=L"r_p = %$(rp), r_c = %$(rc)," * "\n" * L" K = %$(K),c \approx %$(round(c;digits=3)),\alpha = %$(α)")


pl_pp = plot(pp, pp_orig; layout=(1, 2), size=(800, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm,dpi=300)

savefig(pl_pp,figs_dir*"orig_vs_modelled_pp_c_$c"*".png")

