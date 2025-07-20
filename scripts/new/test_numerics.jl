using Revise
using Plots,LaTeXStrings
using AttentionDynamics
#using DifferentialEquations
using OrdinaryDiffEq
using DynamicalSystems

#N = 10
#g = 45.0
#w = 0.5

N = 5
g = 20.0
w = 0.4
p = [g,w,N]

dts = [0.1,0.05,0.001,0.005,0.0001]
algs = [RK4(),Tsit5(),Rodas5(),Vern9()] 
adaptive_settings = [true,false]

T = 100
Ttr = 0
u0 = rand(2N)
colors = distinguishable_colors(length(algs))
alg_strings = ["RK4","Tsit5","Rodas5","Vern9"] 
adaptive_settings_strings = ["adaptive","fixed"]


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
            traj,t = trajectory(ds,T;Ttr=Ttr)
            plot!(pl,t,traj[:,1];lc=colors[i],
                label = alg_strings[i],
                title= adaptive_settings_strings[a] * ", dt=$dt")
        end
    push!(plots,pl)    
    end
end

pl_all = plot(plots...,layout=(length(dts),2),size=(800,1200))