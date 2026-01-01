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
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds
u01 = 0.1 .* rand(2N)
ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

#length of precomputed trajs 
T = 100 #5000
Ttr = 100 #500
Δt = 1e-1

#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

rs = [6:0.1:15;]
#rs4 = [9.0:1.0:12.0;]
#rs4 = [6.0:1.0:9.0;]
#rs4 = [6.0:10.0:16.0;]
#rs4 = [6.0:10.0:40.0;]
#rs4 = [6.0:10.0:36.0;]
rs4 = [6.0:7.0:36.0;]

@show rs4
set_parameters!(ds, p_orig)
rescaled_params.(rs4, rs4, K, c, α, N)

#precompute trajectories 
trajs = precompute_trajectories(ds, rs4, [1, 2], Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_p=r_c = %$(rs4[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"r_p=r_c = %$(rs4[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(length(rs4), 2), size=(600, 800), dpi=300)
#save traj plot 

savedir = "../../data/trajs_precomp_lifetimes_test2/rp=rc/"
mkpath(savedir)
savefig(pl_trajs, savedir*"trajs_x_y_change_rp_rc.png")


#save trajectories into file 
for i in 1:length(trajs)    
    traj = trajs[i]
    writedlm(savedir*"trajs_rp=rc$(rs4[i])_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end

#---------------------------------change rp---------------------------------
# compare with rp \neq rc (mismatch between production and consumption rate)

set_parameters!(ds, p_orig)
rescaled_params.(rs4, rs4, K, c, α, N)

#precompute trajectories 
trajs = precompute_trajectories(ds, rs4, 1, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_c=%$(rc), r_p = %$(rs4[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"r_c=%$(rc), r_p = %$(rs4[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(4, 2), size=(600, 800), dpi=300)
#save traj plot 

savedir = "../../data/trajs_precomp_lifetimes_test2/rp/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_rp.png")


#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rp$(rs4[i])_rc$(rc)_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end

#---------------------------------change rc---------------------------------

#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)
#ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)
set_parameters!(ds, p_orig)
#set_parameter!(ds, 1,6.0)
trajs = precompute_trajectories(ds, rs4[end-1:end], 2, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"r_p = %$(rp), r_c = %$(rs4[i])")
    pl_y = plot(Matrix(traj[1:5000, 2:2:16]), legend=false, title=L"r_p = %$(rp), r_c = %$(rs4[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(4, 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_lifetimes_test2/rc/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_rc.png")


#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rc$(rs4[i])_rp$(rp)_K$(K)_c$(c)_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end



#----------------------------change α parameter-------------------------
#increase α -> decrease memory length 

#------------params---------------

#gs4 = rs4 ./ α
#αs4 = rp ./ gs4

αs4 = [0.01, 0.0085, 0.0025, 0.00175]
rescaled_params.(rp, rp, K, c, αs4, N)

#loop over memory length
#precompute trajectories 
set_parameters!(ds, p_orig)
trajs = precompute_trajectories(ds, αs4, 5, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for traj in trajs
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false)
    pl_y = plot(Matrix(traj[1:8500, 2:2:16]), legend=false)

    push!(plots_trajs, pl_x, pl_y)
end

pl_trajs = plot(plots_trajs..., layout=(4, 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_lifetimes_test2/alpha/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_alpha.png")

#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rp_$(rp)_rc$(rc)_K$(K)_c$(c)_α$(αs4[i])_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end

#----------------------------change c parameter-------------------------
#increase c -> increase competition 


cs = [1.5:0.3:2.5;]
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-1, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)


u01 = 0.1 .* rand(2N)

ds = create_original_ds(u01, deepcopy(p_orig); diffeq=diffeq)

#precompute trajectories 
set_parameters!(ds, p_orig)
trajs = precompute_trajectories(ds, cs, 4, Δt, T, Ttr)

#plot parts of trajectories  
plots_trajs = []

for (i, traj) in enumerate(trajs)
    pl_x = plot(Matrix(traj[1:500, 1:2:2N]), legend=false, title=L"c = %$(cs[i])")
    pl_y = plot(Matrix(traj[1:500, 2:2:16]), legend=false, title=L"c = %$(cs[i])")

    push!(plots_trajs, pl_x, pl_y)
end
pl_trajs = plot(plots_trajs..., layout=(4, 2), size=(600, 800), dpi=300)

savedir = "../../data/trajs_precomp_lifetimes_test2/c/"
mkpath(savedir)
savefig(pl_trajs, savedir * "trajs_x_y_change_c.png")


#save trajectories into file 
for i in 1:length(trajs)
    traj = trajs[i]
    writedlm(savedir * "trajs_rp_$(rp)_rc$(rc)_K$(K)_c$(cs[i])_α$(α)_N$(N)_T$(T)_Ttr$(Ttr)_deltat$(Δt).txt", Matrix(traj))
end