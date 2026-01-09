using Plots, LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
cd(@__DIR__)

include("chart_functions.jl")

k = 10 #top k
Δt_chart = 1 #charting time

#READ precomputed lifetimes
data_dir_lt_root = "../../data/lifetimes_precomp_test7_long/"
data_dir_lt = data_dir_lt_root * "rp=rc/"
files_lt_all = readdir(data_dir_lt)
files_lt = filter(f -> occursin("Δt1.", f) && occursin("k_10_", f), readdir(data_dir_lt))


#READ precomputed trajectories
data_dir_traj_root = "../../data/trajs_precomp_test7_long/"
data_dir_traj = data_dir_traj_root*"rp=rc/"
files_trajs = readdir(data_dir_traj)
traj_files = filter(f -> occursin(".txt", f), files_trajs)[2:end]
#files_lt = filter(f -> occursin("Δt1.", f) && occursin("k_10_", f), readdir(data_dir_lt))
rs = vec(readdlm(data_dir_traj * files_trajs[1]))


trajs = Matrix{Float64}[]

for (i,f) in enumerate(traj_files)
    println("Reading file " * f * "...")
    traj = readdlm(data_dir_traj * f, Float64)
    push!(trajs, traj)
end

#lifetimes
lifetimes_x = [] 
for (i,f) in enumerate(files_lt)
    lts = readdlm(data_dir_lt*f)
    push!(lifetimes_x, lts)
end
mean_τs = mean.(lifetimes_x)

#inequality
mean_ineqs = get_inequality(trajs, k, Δt_chart)
mean_ineq = [ineq[1] for ineq in mean_ineqs ]

#---------------------------------lyapunov_exp----------------------------
N = 300

#original params
#rp=3.0;rc = 3.0;K = 1.0;α = 5*1e-3;c=1.8
rc = 23.0;
rp = 23.0 #9.0;
K = 1.0;
c = 2.4;
α = 0.005;
w, g = rescaled_params(rp, rc, K, c, α, N)

p_orig = [rp, rc, K, c, α]

#critical params from the rescaled 
wc1 = (N - 1) / (rp / α)
wc2 = (N - 1) / N


#------------------------------solver settings----------------------------------
dtmax = 1e-1 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=5*1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds
u01 = 0.1 .* rand(2N)
ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

#lyap calc settings 
#lyapunov calculation is not affected by sticking?
#perturbations kick the system out of sticked states?
#check this hypothesis
T = 2000 #5000
Ttr = 1000 #500
Δt = 1e-1

set_parameters!(ds, p_orig)

λs = []
for (i,r) in enumerate(rs)
    @show r
    set_parameter!(ds, 1, r)
    set_parameter!(ds, 2, r)

    λ =  lyapunov(ds,T;Ttr=Ttr)
    push!(λs, λ)
end

#-----------------test plot---------------------

#plot critical curve
ws = [0.01:0.01:0.99;]

#when rp=rc
w_g_rs = rescaled_params.(rs, rs, K, c, α, N)
w_rs = [w_g_rs[i][1] for i in 1:length(w_g_rs)]
g_rs = [w_g_rs[i][2] for i in 1:length(w_g_rs)]

#mean lifetimes 
pl_mean_lt = plot(ws, 1 ./ ws, lc=:gray10, 
    lw=1, ylabel=L"g/(N-1)", xlabel=L"w",
    guidefontsize=15, label=L"w^H", 
    title=L"\langle \tau \rangle", 
    xlims=(0.0, 0.5), ylims=(0., 25),
    legendfontsize=10)

plot!(pl_mean_lt, w_rs, g_rs ./ (N - 1),
    zcolor=mean_τs, m=(:default, 10.0, Plots.stroke(1, :gray10)),
    lw=1, #arrow=:arrow,
    label=L"r_p=r_c \in [%$(rs[1]),%$(rs[end])]", ms=4,
    markershape=:circle)

#mean ineq 
pl_ineq = plot(ws, 1 ./ ws, lc=:gray10, lw=1, 
    ylabel="", xlabel=L"w",
    guidefontsize=15, label=L"w^H", 
    title="Mean inequality",
    xlims=(0.0, 0.5), ylims=(0., 25),
    yformatter=:none,
    legendfontsize=10)

plot!(pl_ineq, w_rs, g_rs ./ (N - 1),
    zcolor=mean_ineq, m=(:default, 10.0, Plots.stroke(1, :gray10)),
    lw=1, #arrow=:arrow,
    label=L"r_p=r_c \in [%$(rs[1]),%$(rs[end])]", ms=4,
    xlims=(0.0, 0.5), ylims=(0., 25), markershape=:circle)


pl_lyap = plot(ws, 1 ./ ws, lc=:gray10, lw=1, 
    ylabel="", xlabel=L"w",
    guidefontsize=15, label=L"w^H", 
    title=L"\lambda", 
    xlims=(0.0, 0.5), ylims=(0., 25),
    yformatter=:none,
    legendfontsize=10)

#lyapunov exponent 
plot!(pl_lyap, w_rs, g_rs ./ (N - 1),
    zcolor=λs, m=(:default, 10.0, Plots.stroke(1, :gray10)),
    lw=1, #arrow=:arrow,
    label=L"r_p=r_c \in [%$(rs[1]),%$(rs[end])]",
    ms=4,xlims=(0.0, 0.5), ylims=(0., 25), 
    markershape=:circle)

pl = plot(pl_mean_lt,pl_ineq,pl_lyap,
    layout=(1,3),dpi=300,size=(1800,500),
    xlims=(0.1,0.5),
    plot_title=L"N=%$(N), K=%$(K), c=%$(c)")

savedir = "../../paper_figs/"
savefig(pl,savedir*"phase_diag_measures_rp=rc_test7_long.png")