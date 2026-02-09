using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)


function precompute_trajectories(ds, ps, p_idx, Δt=1e-1, T_traj=5000, Ttr=1000)

    #precomputed trajs
    trajs = []
    for p in ps

        if p_idx isa Vector
            for i in p_idx
                set_parameter!(ds, i, p)
            end
        else
            set_parameter!(ds, p_idx, p)
        end

        #solve systems to get trajectory
        @time traj, times = trajectory(ds, T_traj; Ttr=Ttr, Δt=Δt)
        push!(trajs, traj)
    end
    return trajs
end

#-----------------------------------time params---------------------------------------

#length of precomputed trajs 
T = 100000 #10000 #100 #5000
Ttr = 500 #500
Δt = 1e-1 #5*1e-2

#-----------------------------------setup system--------------------------------------

diffeq2 = (alg=Tsit5(), dt=1e-2, dtmax=1e-1, adaptive=true, maxiters=1e15)

N = 10
g = 18.0
w = 0.7
u0 = rand(2N) * 0.1

ds = create_ds(u0, [g, w, N]; diffeq=diffeq2) #no noise

#parameters to change 
gs = [21.11111111111111, 15.200000000000001]
ws = [0.5612009237875288, 0.6398104265402844]

#plot parts of trajectories  
plots_trajs = []

#precompute trajectories
savedir = "../../data/trajs_precomp_rescaled_alpha/"
writedlm(savedir*"g_params.txt", gs)
writedlm(savedir*"w_params.txt", ws)

mkpath(savedir)
plots_trajs = []

for i in 1:length(ws)
    @show i
    #set parameters
    w = ws[i]
    g = gs[i]

    set_parameter!(ds, 1, g)
    set_parameter!(ds, 2, w)

    #solve systems to get trajectory
    @time traj, times = trajectory(ds, T; Ttr=Ttr, Δt=Δt)

    pl = plot(Matrix(traj[1:500, 1:2:2N]), title=L"w=%$(w),g=%$(g)",legend=false)
    push!(plots_trajs,pl)

    xs = vcat(Matrix(traj)[:,1:2:end]...)
    write(savedir * "trajs_g_$(g)_w$(w)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat", xs)

end

pl_trajs = plot(plots_trajs..., layout=(2, 1), size=(1000, 300), dpi=300)
savefig(pl_trajs, savedir * "trajs_x_change_alpha.png")
