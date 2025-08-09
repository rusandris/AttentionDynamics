using Revise
using Plots,LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
cd(@__DIR__)
#inspectdr()
gr()

function plot_projection(traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf))
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    pl = plot(x1 .- x2, y1 .- y2,
        lw=1, la=1.0, ylims=ylims, xlims=xlims)

end

function plot_projection!(pl,traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf),label="")
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    plot!(pl,x1 .- x2, y1 .- y2,
        lw=1, la=1.0, ylims=ylims, 
        xlims=xlims,label=label)

end


#check if partial limit cycles really exist
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.3 #w = 0.25

p = [g,w,N]

#search a u0 before transition to limit cycle
T = 20
Ttr = 2000



#plot usual and stiff solve on top of each other many times

plots = []
for _ in 1:10
    u0 = rand(2N) .* 50
    #default diffeq options
    diffeq0 = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
    ds0 = create_ds(u0,p;diffeq=diffeq0)
    traj0,t0 = trajectory(ds0,T,u0;Ttr=Ttr,Δt = 0.001)
    M0 = Matrix(traj0)
    pl = plot_projection(traj0;xlims=(-30,30),ylims=(-4,4))
    #plot(t,M[:,5],xlims=(0,100))


    adaptive_settings_strings = ["adaptive","fixed"]
    Δt = 1e-2

    #----------------solver options------------------ 

    dts = [1e-4] #, 1e-4] #,1e-5] #,1e-6]
    algs = [Rodas5()] #[Rosenbrock23(),RK4(),Tsit5(),Rodas5(),Vern9()] 
    alg_strings = ["Rodas5"] #["Rosenbrock23","RK4","Tsit5","Rodas5","Vern9"] 

    adaptive_settings = [true,false]

    #calc trajectories with different solvers and dt
    trajs_all = []
    for dt in dts
        @show dt
        trajs = []
        for (j,alg) in enumerate(algs)
            @show alg_strings[j]
            diffeq = (alg = alg,dt = dt,dtmax=dt,adaptive=true,maxiters=1e15)
            ds = create_ds(u0,p;diffeq=diffeq)
            traj,t = trajectory(ds,T,u0;Ttr=Ttr,Δt = 0.01)
            push!(trajs,traj)
        end
        push!(trajs_all,trajs)
    end

    #plot on top of eachother
    colors = distinguishable_colors(length(algs))
    for traj in trajs_all[1] 
        plot_projection!(pl,traj;xlims=(-30,30),ylims=(-4,4),label=alg_strings[1])
    end
    push!(plots,pl)
end

pl2 = deepcopy(pl)
for traj in trajs_all[2] 
    plot_projection!(pl2,traj;xlims=(-30,30),ylims=(-4,4))
end
pl2

#=
plots = []
for (a,setting) in enumerate(adaptive_settings)
    @show setting
    for dt in dts
        @show dt
        pl = plot()
        for (i,alg) in enumerate(algs)
            @show alg_strings[i]
            
            diffeq = (alg=alg,dt=dt,dtmax=dt,adaptive=setting)
            ds = create_ds(u0,p;diffeq)
            traj,t = trajectory(ds,T;Ttr=Ttr,Δt = Δt)
            plot!(pl,traj[:,1],traj[:,2];lc=colors[i],
                label = alg_strings[i],
                title= adaptive_settings_strings[a] * ", dt=$dt"*" T=$T"*" Ttr=$Ttr")
        end
    push!(plots,pl)    
    end
end

pl_all = plot(plots...,layout=(length(dts),2),size=(800,1200))
savefig(pl_all,"../../new_figs/solver_tests_N_$(N)_w_$(w)_g_$(g)"*".png")
=#

#=
#----------------calculate lyapunov----------------------
Ttr = 2000
T = 2000
for (a,setting) in enumerate(adaptive_settings)
    @show setting
    for dt in dts
        @show dt
        pl = plot()
        for (i,alg) in enumerate(algs)
            @show alg_strings[i]
            
            diffeq = (alg=alg,dt=dt,dtmax=dt,adaptive=setting)
            ds = create_ds(u0,p;diffeq)
            λ = lyapunov(ds,T;Ttr=Ttr)
            @show (setting,dt,alg_strings[i], λ)

        end

    end

end
=#
