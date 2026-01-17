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


#setup system 

N = 300

#original params
#rp=3.0;rc = 3.0;K = 1.0;α = 5*1e-3;c=1.8
r_const = 10.0
rc = r_const
rp = r_const #9.0;
K = 1.0;
c = 2.4;
α = 0.005;
w, g = rescaled_params(rp, rc, K, c, α, N)

p_orig = [rp, rc, K, c, α]
#critical params from the rescaled 
wc1 = (N - 1) / (rp / α)
wc2 = (N - 1) / N


#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds

#u01 = 0.1 .* rand(2N)
#u01 = BigFloat.(u0;precision=56) #turn to BigFloat (higher precision)
x0 = rand(N) * 0.001
u0 = zeros(2N)
u0[1:2:end] .= x0

ds = create_original_ds(u0, deepcopy(p_orig); diffeq=diffeq)

σ = 1e-12
p_orig_noise = [rp, rc, K, c, α,σ]
ds_noise = create_original_ds_noisy(u0, p_orig_noise; diffeq=diffeq) #noise 

#length of precomputed trajs 
T = 100 #10000 #100 #5000
Ttr = 1000 #500
Δt = 1e-1

#----------------------------change r parameter-------------------------
rs = [6.0:4.0:24.0;]

@show rs
set_parameters!(ds, p_orig)
rescaled_params.(rs, rs, K, c, α, N)


#precompute trajectories 
trajs = precompute_trajectories(ds_noise, rs, [1, 2], Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:50, 1:2:2N]), legend=false, title=L"r_p=r_c = %$(rs[i])")
    pl_y = plot(Matrix(traj[1:50, 2:2:16]), legend=false, title=L"r_p=r_c = %$(rs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(rs), 2), size=(600, 800), dpi=300)
#save traj plot 