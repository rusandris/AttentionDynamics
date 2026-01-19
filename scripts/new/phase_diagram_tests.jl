using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
using LinearAlgebra,StatsBase
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")


function trajpair_distance(ds, n_ens, nr_steps; d0=1e-12, Ttr_init=3000, T_init=100)

    traj_u0s, = trajectory(ds, T_init; Ttr=Ttr_init) #generate init conds 
    u0s = rand(traj_u0s, nr_ens)
    pds = ParallelDynamicalSystem(ds, [u0s[1], u0s[1]]) #dummy init

    #calc d(t) avgs
    dist_avg = zeros(nsteps)
    for (i, u0) in enumerate(u0s)
        reinit!(pds, [u0, u0 .+ rand(2N) .* d0 / sqrt(2N)])

        d = zeros(nsteps)
        for i in 1:nsteps
            step!(pds)
            states = current_states(pds)
            dist_avg[i] += norm(states[1] .- states[2])
        end
    end
    dist_avg ./= nr_ens
    return dist_avg
end

#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients

#read w and g params
#ws_saved = readdlm("../../plotting/data/N=10.0_wvalues.txt")
#gs_saved = readdlm("../../plotting/data/N=10.0_gvalues.txt")
#gs = gs_saved / (10-1)

#read param values
#data_dir = "./phasediag_data_N10/"
data_dir = "../../data/phasediag_data_N10_zoomed/"
#mkpath(data_dir)
figs_dir = "../../paper_figs/"
#mkpath(figs_dir)
saved_files = readdir(data_dir)
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[1:end,1] #cut some off the end 
gs = params[1:end,2]


#=
ws = range(0.5,0.85;length=50)
g_rel = 5
gs = range(1.2,3.0;length=50) .* (N-1)
=#
#=
#try original equations
u0 = generate_initial_condition(N)
rc=9.0;rp=9.0;K=1.0;c=2.4;α=0.005
p = [rc,rp,K,c,α]
ds = create_original_ds(u0,p)
=#

#critical params
wc1 = (N-1)/gs[1]
wc2 = (N-1)/N


#---------------------time params----------------------
T_fp = 100
Ttr_fp = 100
Ttr_init = 3000
Δt_fp = 0.1
fp_settings = [T_fp,Ttr_fp,Δt_fp]
Ttr = 2000
Ttrs = [Ttr]
T = 5000
lyap_vals = []
d0=1e-12
#calc lyapunov exponents
@time for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds_noise,gs,ws; T=T, Ttr=Ttr, Δt=1.0, fp_settings=fp_settings, 
        d0=d0, ut=1e-3, lt=1e-14,save=false,realtime_plotting=false,parallelize=false)

    #@time Λ = lyapunov_diagram(ds,gs,ws; T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
    #d0=1e-9,save=false,realtime_plotting=false,parallelize=false)
    push!(lyap_vals,Λ)
end

savedir = "../../data/phasediag_data_N10_noise_$(σ)_d0_$(d0)_ut1e-3/"
mkpath(savedir)
writedlm(savedir * "phasediag_params.txt", hcat(ws, gs))
writedlm(savedir * "lyap_values.txt", lyap_vals)


#--------------------------check lyap-----------------------
#average trajectory distance evolution 
#no noise vs noise case 

#no noise
n_pairs = 3
nsteps = 1500
N = 10
g = 18.0
w = 0.7
ts = [1:nsteps+1;] ./ 1e2

#diffeq = (alg=Tsit5(), dt=1e-2, adaptive=false, maxiters=1e15)
diffeq = (alg=RK4(), dt=1e-2, adaptive=false, maxiters=1e15)
u0 = rand(2N) * 0.1
#----------------------systems----------------------
ds = create_ds(u0, [g, w, N]; diffeq=diffeq) #no noise
#ds = Systems.lorenz()
#noise
σ = 1e-6
ds_noise = create_ds_noisy(u0, [g, w, N, σ]; diffeq=diffeq) #noise 

dim = dimension(ds)


#initial conditions
traj_u0s, = trajectory(ds, 100; Ttr=Ttr_init) #generate init conds 
u0s = rand(traj_u0s, n_pairs)

#TODO: make a summary plot for initial distance scaling
#chaotic and regular flow with two initial distances 
#add regular flow with noise 

d0s = [1e-16,1e-12,1e-5]
dmins = Float64[]
dmaxs = Float64[]
plots = []
σs = [1e-12, 1e-6]
ws = [0.6,0.7]
dists_all_sys = []

for ds in [ds,ds_noise]
    dists_params = []
    for w in ws 
        dists_d0s = []
        @show w
        set_parameter!(ds,2,w)
        pds = ParallelDynamicalSystem(ds, [u0, u0]) #dummy init
        for d0 in d0s
            dists_all_pairs = Vector{Float64}[]
            @show d0
            dist_avg = zeros(nsteps)
            for (i,u0) in enumerate(u0s)

                @show i
                init_states = [u0, u0 .+ rand(dim) .* d0 / sqrt(dim)]
                reinit!(pds, init_states)
                d_init = norm(init_states[1] .- init_states[2])

                dists = zeros(nsteps+1)
                dists[1] = d_init
                for i in 1:nsteps
                    step!(pds)
                    states = current_states(pds)
                    d = norm(states[1] .- states[2])

                    #save every dist for every pair 
                    dists[i+1] = d

                end
                push!(dists_all_pairs,dists)
            end
            dist_avg ./= n_pairs

            push!(dists_d0s, dists_all_pairs)

            #d_pair_maxs = [maximum(d[end-1000:end]) for d in dists_all_pairs]
            #dmax,dmax_idx = findmax(d_pair_maxs)
            #dmin, dmin_idx = findmin(d_pair_maxs)
            #@show dmax_idx, dmin_idx
            #push!(dmins,dmin)
            #push!(dmaxs, dmax)
        end
        push!(dists_params, dists_d0s)
    end
    push!(dists_all_sys, dists_params)
end

#test plot
pl = plot(legend=false, lw=2, yticks = 10.0 .^ [-18:2:2;],
    yscale=:log10, ylabel=L"d_{12}(t)", xlabel=L"t",
    legendfontsize=12, ylims=(1e-19, 100),
    title=L"w=%$(ws[1]), \sigma=0")
    #label=L"N = %$N, g=%$g, w=%$w")

plots_orig = [pl, deepcopy(pl)]
plots_noise = [deepcopy(pl), deepcopy(pl)]

plot!(plots_orig[2], title=L"w=%$(ws[2]), \sigma=0",ylabel="",yformatter=:none)
plot!(plots_noise[1], title=L"w=%$(ws[1]), \sigma=%$σ")
plot!(plots_noise[2], title=L"w=%$(ws[2]), \sigma=%$σ", ylabel="", yformatter=:none)
plots = [plots_orig,plots_noise]

colors = [:gray10, :gray30, :gray50]

for (s,dists_sys) in enumerate(dists_all_sys)
    pls = plots[s]
    for (i,w) in enumerate(ws)

        pl = pls[i]

        dists_p = dists_sys[i]
        for (j,d0) in enumerate(d0s)
            dists_d0 = dists_p[j]
            for (k,u0) in enumerate(u0s)
                dists = dists_d0[k]

                plot!(pl,ts, dists, lc=colors[k],lw=1,
                    label=L"N = %$N, g=%$g, w=%$w")

            end 
        end 
    end
end
pl2 = plot(plots[1]...,plots[2]...,layout=(2,2),size=(600,400),dpi=300)
savefig(pl2,"../../paper_figs/initital_distance_scaling_test.png")
#TODO: order parameter (λ) versus noise amplitude
#including \sigma = 0 and two (one long and a short transient time)

σs = [10.0 .^ [-14:1:-1;]...]
lambdas = Float64[]
d_long = Float64[]
nr_ens = 1000
nsteps = 2000
interval = 1:300
ts = [1:nsteps;] ./1e2
pl = plot()
for σ in σs
    @show σ
    set_parameter!(ds_noise,4,σ)
    #λ = lyapunov(ds_noise, 1000; Ttr=0, u0=u0s[1], d0_lower=1e-14, d0_upper=1e-8, d0=1e-12, Δt=1.0)

    d = trajpair_distance(ds_noise, nr_ens, nsteps; Ttr_init=3000, T_init=3000)
    long_term_mean = mean(d[end-500:end])
    plot!(ts,d)
    #itsct, λ = linreg(ts[interval], log.(d[interval]))

    push!(d_long,long_term_mean)
end

pl = plot(σs,d_long,ylabel=L"d(t>>)",
    xlabel=L"\sigma",xscale=:log10,
    yscale=:log10, xticks=σs ,
    legendfontsize=12,
    guidefontsize=22,
    label=L"N=%$(N),w=%$(w),g=%$(g),\sigma=%$(σ)",
    legend=:bottomright,
    st=:scatter,ms=3,mc=:gray10,dpi=300)

savefig(pl, "../../paper_figs/long_term_average_distance_vs_noise_scaling_test.png")

#=
#check slopes 
ts = [1:nsteps;] ./ 1e2
interval = 1:300

ts_int = ts[interval]
#noise
itsct, slope = linreg(ts_int,log10.(dist_avg_noise[interval]))
λ = slope * log(10)
plot!(ts[interval], 10 .^ (itsct .+ slope .* ts_int), 
    lc=:gray30, label=L"\lambda = %$λ")

#noise
interval = 300:3000
ts_int = ts[interval]
itsct, slope = linreg(ts_int, log10.(dist_avg_noise[interval]))
λ = slope * log(10)
plot!(ts[interval], 10 .^ (itsct .+ slope .* ts_int),
    lc=:gray30, label=L"\lambda = %$λ")

#no noise
interval = 1:300
ts_int = ts[interval]
itsct, slope = linreg(ts_int, log10.(dist_avg[interval]))
λ = slope * log(10)
plot!(ts_int, 10 .^ (itsct .+ slope .* ts_int),
    lc=:orange, label=L"\lambda = %$λ")

#no noise
interval = 2000:6000
ts_int = ts[interval]
itsct, slope = linreg(ts_int, log10.(dist_avg[interval]))
λ = slope * log(10)
plot!(ts_int, 10 .^ (itsct .+ slope .* ts_int),
    lc=:orange, label=L"\lambda = %$λ")


#--------------read lyaps--------------------
data_dir = "../../data/phasediag_data_N10_noise_1e-4/"
lyap_vals = reshape(readdlm(data_dir*"lyap_values.txt"),(50,50))
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[:,1]
gs = params[:,2]
N=10

pl_lyaps = plot_phase_diagram(lyap_vals[1],ws,gs,N;HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}",clims=(-0.5,0.5),
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large)
#pl2 = plot_phase_diagram(lyap_vals[2],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
#pl3 = plot_phase_diagram(lyap_vals[3],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
#savefig(pl1,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])_noise_1e-2"*".png")
savefig(pl_lyaps,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])"*".png")

#savefig(pl2,"../../N5_phasediagram_Ttr_max_$(Ttrs[2])"*".png")
#savefig(pl3,"../../N5_phasediagram_Ttr_max_$(Ttrs[3])"*".png")


plot(ws,lyap_vals,ylims=(-0.2,Inf),labels = ["Ttrs=$(Ttrs[1])" "Ttrs=$(Ttrs[2])" "Ttrs=$(Ttrs[3])"])
=#