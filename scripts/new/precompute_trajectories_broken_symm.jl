using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)


N = 300
K = 1.0
c = 2.4
r = 10.0
α0 = 0.005
#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)

T = 10
Ttr = 0
Δt_traj = 0.1
nbins = 20

#create classic ds
x0 = rand(N) * 0.01
u0 = zeros(2N)
u0[1:2:end] .= x0

#------------------try asymmetric coupling---------------------

#setup interaction matrix
C = randn(N, N) .* 0.5 .+ 2.0
for i in 1:N
    C[i, i] = 0.0
    null_idxs = rand(1:N, rand(50:250))
    C[i, null_idxs] .= 0.0
    #anti_idxs = rand(1:N, rand(1:50))
    #C[i, anti_idxs] .= C[i, anti_idxs] .* (-1.0)
end

ds_as = create_original_asymmetric_interaction(u0, [r, r, K, C, α0]; diffeq=diffeq)

rs = [6.0,8.0,10.0]
p_sets_r = [[rs[i], rs[i]] for i in 1:length(rs)]
nr_ics = 10
ic_counter = 0 
savedir_traj_root = "../../data/trajs_precomp_asymmetric/"
savedir_traj = savedir_traj_root * "asymmetric_coupling/"
mkpath(savedir_traj)

#test precompute traj
precompute_trajectories(ds_as, nr_ics, p_sets_r, [1,2]; Δt=Δt_traj,
    T_traj=T, Ttr=Ttr, savedir=savedir_traj, p_names=["r"],
    ic_counter=ic_counter,
    filename="asymmetric_coupling_N$(N)_K$(K)_alpha$(α0)_T$(T)_Ttr$(Ttr)_deltat$(Δt_traj).dat")


#------------------try asymmetric coupling with varied delay---------------------

#setup interaction matrix
C = randn(N, N) .* 0.5 .+ 2.0
for i in 1:N
    C[i, i] = 0.0
    null_idxs = rand(1:N, rand(50:250))
    C[i, null_idxs] .= 0.0
    #anti_idxs = rand(1:N, rand(1:50))
    #C[i, anti_idxs] .= C[i, anti_idxs] .* (-1.0)
end
αs = rand(N) .* 0.01 .+ 1e-3

ds_asd = create_original_asymmetric_interaction_varied_delay(u0, [r, r, K, C, αs]; diffeq=diffeq)

rs = [6.0, 8.0, 10.0]
p_sets_r = [[rs[i], rs[i]] for i in 1:length(rs)]
nr_ics = 10
savedir_traj_root = "../../data/trajs_precomp_asymmetric/"
savedir_traj = savedir_traj_root * "asymmetric_coupling_varied_delay/"
mkpath(savedir_traj)

#test precompute traj
precompute_trajectories(ds_as, nr_ics, p_sets_r, [1, 2]; Δt=Δt_traj,
    T_traj=T, Ttr=Ttr, savedir=savedir_traj, p_names=["r"],
    ic_counter=ic_counter,
    filename="asymmetric_coupling_varied_delay_N$(N)_K$(K)_T$(T)_Ttr$(Ttr)_deltat$(Δt_traj).dat")
