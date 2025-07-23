using Revise
using Plots,LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems
cd(@__DIR__)
#N = 10
#g = 45.0
#w = 0.5

N = 5
g = 20.0
w = 0.25
p = [g,w,N]

dts = [0.1,0.05,0.001,0.005,0.0001]
algs = [RK4(),Tsit5(),Rodas5(),Vern9()] 
adaptive_settings = [true,false]

T = 200
Ttr = 0
u0 = rand(2N)
colors = distinguishable_colors(length(algs))
alg_strings = ["RK4","Tsit5","Rodas5","Vern9"] 
adaptive_settings_strings = ["adaptive","fixed"]
Δt = 1e-3


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
