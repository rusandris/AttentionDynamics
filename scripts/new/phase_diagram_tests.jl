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

diffeq = (alg = Tsit5(),dt=1e-2,adaptive=false,maxiters=1e15)
diffeq = (alg = RK4(),dt=1e-2,adaptive=false,maxiters=1e15)

N = 10
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

#-------------------------create ds---------------------------------
#generate standard init cond
#u0 = rand(2N)
u0 = generate_initial_condition(N)
ds = create_ds(u0,[5,0.5,N];diffeq=diffeq) #no noise

σ = 1e-6
p_noisy = [5, 0.5, N, σ]
ds_noise = create_ds_noisy(u0,p_noisy;diffeq=diffeq) #noise 
#figure out how to calc lyap exp for noisy ds 
#d0_lower,d0_upper should be adjsted 



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
n_trajs = 500
d0 = 1e-4
nsteps = 65000
g = 18.0
w = 0.7
ts = [1:nsteps;] ./ 1e2

#----------------------systems----------------------
ds = create_ds(u0, [g, w, N]; diffeq=diffeq) #no noise
#ds = Systems.lorenz()
#noise
σ = 1e-3
ds_noise = create_ds_noisy(u0, [g, w, N, σ]; diffeq=diffeq) #noise 

dim = dimension(ds)
pds = ParallelDynamicalSystem(ds_noise, [u0, u0]) #dummy init

traj_u0s, = trajectory(ds, 100; Ttr=Ttr_init) #generate init conds 
u0s = rand(traj_u0s, n_trajs)
#calc d(t) avgs
d0s = [1e-12,1e-10,1e-8,1e-6,1e-4,1e-2]
dmins = Float64[]
dmaxs = Float64[]
plots = []
σs = [1e-12, 1e-10, 1e-8, 1e-6]
d0 = 1e-12
for σ in σs
    @show σ

    set_parameter!(pds,4,σ)
    dist_avg = zeros(nsteps)
    dists_all_pairs = Vector{Float64}[]

    for (i,u0) in enumerate(u0s)

        #@show i
        reinit!(pds, [u0, u0 .+ rand(dim) .* d0 / sqrt(dim)])
        dists = zeros(nsteps)
        d = zeros(nsteps)
        for i in 1:nsteps
            step!(pds)
            states = current_states(pds)
            d = norm(states[1] .- states[2])

            #save every dist for every pair 
            dists[i] = d

            #add to avg  
            dist_avg[i] += d
        end
        push!(dists_all_pairs,dists)
    end
    dist_avg ./= n_trajs

    d_pair_maxs = [maximum(d[end-1000:end]) for d in dists_all_pairs]

    dmax,dmax_idx = findmax(d_pair_maxs)
    dmin, dmin_idx = findmin(d_pair_maxs)

    @show dmax_idx, dmin_idx

    #test plot 
    pl = plot(dists_all_pairs[dmax_idx][end-2000:end],
        lc=:red,yscale=:log10,
        label=L"d_{max}(t)")

    plot!(dists_all_pairs[dmin_idx][end-2000:end],
        lc=:blue,
        label=L"d_{min}(t)")

    push!(plots,pl)

    push!(dmins,dmin)
    push!(dmaxs, dmax)
end

plot!(ts,dist_avg, lc=:lightblue, yscale=:log10,
    label=L"\sigma = 0, d_0 = %$d0",
    title = L"g=%$g,w=%$w",
    legend=:bottomright,lw=2,ticks=10,
    legendfontsize=12)


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