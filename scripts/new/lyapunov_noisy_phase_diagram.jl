using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots, LaTeXStrings
using LinearAlgebra
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")

function lyapunov_from_slope(ds, n_ens, nr_steps;d0=1e-12, Ttr_init=3000, T_init=100)

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

#read param values
#data_dir = "./phasediag_data_N10/"
data_dir = "../../data/phasediag_data_N10_zoomed/"
#mkpath(data_dir)
figs_dir = "../../paper_figs/"
#mkpath(figs_dir)
saved_files = readdir(data_dir)
params = readdlm(data_dir * "phasediag_params.txt")
ws = params[1:end, 1] #cut some off the end 
gs = params[1:end, 2]

diffeq = (alg=Tsit5(), dt=1e-2, adaptive=false, maxiters=1e15)
diffeq = (alg=RK4(), dt=1e-2, adaptive=false, maxiters=1e15)

N = 10

#TODO: automate slopefitting to get λ from noisy system

σ = 1e-6
p_noisy = [5, 0.5, N, σ]
ds_noise = create_ds_noisy(u0, p_noisy; diffeq=diffeq) #noise 

interval = 1:300
nr_ens = 100
d0 = 1e-12
nsteps = 15000
ts = [1:nsteps;] ./ 1e2
ts_int = ts[interval]

for (i, g) in enumerate(g)
    @show g
    wc1 = (N - 1) / g

    set_parameter!(ds_noise, 1, g)
    for (j,w) in enumerate(ws)
        w < wc1 && continue
        @show w
        set_parameter!(ds_noise, 2, w)
        d = lyapunov_from_slope(ds_noise, nr_ens, nsteps; Ttr_init=3000, T_init=100)
        
        itsct, λ = linreg(ts_int, log.(d[interval]))
        @show λ
    end
end

plot!(pl, ts, dist_avg_noise, lc=:red, yscale=:log10,
    label=L"\sigma = %$σ, d_0 = %$d0",
    legend=:bottomright, lw=2,
    legendfontsize=12, dpi=300)