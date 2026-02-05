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
r_const = 9.0
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
u02 = zeros(2N)
u02[1:2:end] .= x0

ds = create_original_ds(u02, deepcopy(p_orig); diffeq=diffeq)

#length of precomputed trajs 
T = 10000 #10000 #100 #5000
Ttr = 500 #500
Δt = 1e-1 #5*1e-2

#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

#rs = [6:0.1:15;]
#rs = [6.0:1.0:12.0;]
#rs = [6.0:3.0:16.0;]
rs = [9.0:1.0:12.0;]

@show rs
set_parameters!(ds, p_orig)
rescaled_params.(rs, rs, K, c, α, N)


#precompute trajectories 
trajs = precompute_trajectories(ds, rs, [1, 2], Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_p=r_c = %$(rs[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"r_p=r_c = %$(rs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(rs), 2), size=(600, 800), dpi=300)
#save traj plot 

savedir = "../../data/trajs_precomp_test8_origpar/rp=rc/"
mkpath(savedir)
savefig(pl_trajs, savedir*"trajs_x_y_change_rp_rc.png")
writedlm(savedir*"rp_rc.txt", rs)

#save trajectories into file 
for i in 1:length(trajs)    
    traj = trajs[i]
    writedlm(savedir*"trajs_rp=rc$(rs[i])_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end


#---------------------------------change rp---------------------------------
# compare with rp \neq rc (mismatch between production and consumption rate)
rs_rp = rs
set_parameters!(ds, p_orig)
rescaled_params.(rs_rp, r_const, K, c, α, N)

#precompute trajectories 
trajs = precompute_trajectories(ds, rs_rp, 1, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_c=%$(rc), r_p = %$(rs[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"r_c=%$(rc), r_p = %$(rs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(rs), 2), size=(600, 800), dpi=300)
#save traj plot 

savedir = "../../data/trajs_precomp_test8_origpar/rp/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_rp.png")
writedlm(savedir * "rp_rc=$(r_const).txt", rs_rp)

#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rp$(rs[i])_rc$(rc)_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end

#---------------------------------change rc---------------------------------

#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)
#ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)
r_const = 6.0
set_parameters!(ds, p_orig)
set_parameter!(ds, 1, r_const)

rs_rc = [4.0:0.5:6.0;]


trajs = precompute_trajectories(ds, rs_rc, 2, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_p = %$(rp), r_c = %$(rs[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"r_p = %$(rp), r_c = %$(rs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(rs), 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_test8_origpar/rc/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_rc.png")
writedlm(savedir * "rc_rp=6.txt", rs_rc)

#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rc$(rs[i])_rp$(rp)_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end



#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#------------params---------------

#gs4 = rs ./ α
#αs4 = rp ./ gs4

αs = [0.001, 0.005, 0.008, 0.01]
rescaled_params.(rp, rp, K, c, αs, N)

#loop over memory length
#precompute trajectories 
set_parameters!(ds, p_orig)
trajs = precompute_trajectories(ds, αs, 5, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for traj in trajs
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false)
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false)

    push!(plots_trajs, pl_x, pl_y)
end

pl_trajs = plot(plots_trajs..., layout=(length(αs), 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_lifetimes_alpha/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_alpha.png")

#save trajectories into file 
for i in 1:length(trajs)
    traj = vcat(Matrix(trajs[i])[:,1:2:600]...)
    write(savedir * "trajs_rp_$(rp)_rc$(rc)_K$(K)_c$(c)_α$(αs[i])_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).dat", traj)
end

#----------------------------change c parameter-------------------------
#increase c -> increase competition 

#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#range of c
cs = [1.4:0.2:2.4;]

#other params
K = 1.0;
rp = 10;rc=10;
α = 0.005;

u01 = 0.1 .* rand(2N)
p = [rp, rc, K, c, α]

ds = create_original_ds(u01, deepcopy(p); diffeq=diffeq)

#precompute trajectories 
trajs = precompute_trajectories(ds, cs, 4, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"c = %$(cs[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"c = %$(cs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(cs), 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_test_c/c_low/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_c.png")
writedlm(savedir * "cs.txt", cs)


#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rp_$(rp)_rc$(rc)_K$(K)_c$(cs[i])_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end

#------------------try varied delay system---------------------
αs = rand(1e-6:1e-7:1e-2, 300)
ds_vd = create_original_varied_delay(u02, (rp, rc, K, 4.4, αs); diffeq=diffeq)
traj, times = trajectory(ds_vd, T; Ttr=Ttr, Δt=Δt)
traj_test, times_test = trajectory(ds_vd, 100; Ttr=100, Δt=0.1)

plot(times_test[1:500], Matrix(traj_test[1:500, 1:2:2N]),legend=false)
savedir = "../../data/trajs_precomp_test_varied_delay/1e-6_1e-2/"
mkpath(savedir)
writedlm(savedir * "trajs_rp_$(rp)_rc$(rc)_K$(K)_c$(c)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))


#------------------try varied rate system---------------------
rps = rand(6.0:1e-2:30.0, 300)
ds_vr = create_original_varied_rate(u02, (rps, rc, K, c, α); diffeq=diffeq)
traj, times = trajectory(ds_vr, T; Ttr=Ttr, Δt=Δt)
plot(times[1:5000], Matrix(traj[1:5000, 1:2:2N]), legend=false)
savedir = "../../data/trajs_precomp_test_varied_rate/6_20/"
mkpath(savedir)
writedlm(savedir * "trajs_rp_rc$(rc)_K$(K)_c$(c)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))

#------------------try varied rate and delay system---------------------
rps = rand(6.0:1e-2:30.0, 300)
αs = rand(1e-6:1e-7:1e-1, 300)

ds_vrd = create_original_varied_rate_delay(u02, (rps, rc, K, c, αs); diffeq=diffeq)
traj, times = trajectory(ds_vrd, T ; Ttr=Ttr, Δt=Δt)
traj_test, times_test = trajectory(ds_vrd, 300; Ttr=500, Δt=0.1)
plot(times_test[1:800], Matrix(traj_test[1:800, 1:2:2N]), legend=false)
savedir = "../../data/trajs_precomp_test_varied_rate_delay/"
mkpath(savedir)
writedlm(savedir * "trajs_rc$(rc)_K$(K)_c$(c)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))

#------------------try varied rate and delay system---------------------
#rps = rand(6.0:1e-2:30.0, 300)
rps = fill(9.0,300)
#αs = rand(1e-6:1e-7:1e-1, 300)
αs = fill(α, 300)
cs = [fill(1.1,150)...,fill(4.4,150)...]

ds_vp = create_original_varied_params(u02, (rps, rc, K, cs, αs); diffeq=diffeq)
traj, times = trajectory(ds_vp, T; Ttr=Ttr, Δt=Δt)
traj_test, times_test = trajectory(ds_vp, 230; Ttr=300, Δt=0.1)
plot(times_test, Matrix(traj_test[:, 1:2:2N]), legend=false)
savedir = "../../data/trajs_precomp_test_varied_rate_delay/"
mkpath(savedir)
writedlm(savedir * "trajs_rc$(rc)_K$(K)_c$(c)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
