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

function plot_projection!(pl,traj;xlims=(-Inf,Inf),ylims=(-Inf,Inf),label="",title="")
    x1 = traj[:,1]
    x2 = traj[:,3]
    y1 = traj[:,2]
    y2 = traj[:,4]

    plot!(pl,x1 .- x2, y1 .- y2,
        lw=1, la=1.0, ylims=ylims, 
        xlims=xlims,label=label,title=title)

end


#check if partial limit cycles really exist
N = 10 #N = 5
g = 45.0 #g = 20.0
w = 0.3 #w = 0.25

p = [g,w,N]

#search a u0 before transition to limit cycle
T = 50
Ttr = 1000
nr_traj = 1


#plot usual and stiff solve on top of each other many times

adaptive_settings_strings = ["adaptive", "fixed"]
Δt = 1e-2

#----------------solver options------------------ 

dts = [1e-3, 1e-4] #, 1e-4] #,1e-5] #,1e-6]
algs = [Tsit5(), Rodas5()] #[Rosenbrock23(),RK4(),Tsit5(),Rodas5(),Vern9()] 
alg_strings = ["Tsit5","Rodas5"] #["Rosenbrock23","RK4","Tsit5","Rodas5","Vern9"] 
adaptive_settings = [true, false]
plots = []


#random initial condition
#same for every case
u0 = rand(2N) .* 50
#try adaptive and fixed timestepping
for adapt in adaptive_settings
    plots_adapt = []
    @show adapt
    #calc trajectories with different solvers and dt
    for dt in dts
        @show dt
        #plot it multiple number of times

       
        #Tsit5
        pl0 = plot()
        diffeq0 = (alg=Tsit5(), dt=dt, dtmax=dt, adaptive=adapt, maxiters=1e15)
        ds0 = create_ds(u0, p; diffeq=diffeq0)

        #Rodas5
        pl = plot()
        diffeq = (alg=Rodas5(), dt=dt, dtmax=dt, adaptive=adapt, maxiters=1e15)
        ds = create_ds(u0, p; diffeq=diffeq)

        @time for _ in 1:nr_traj
            #Tsit5
            u0_T5 = current_state(ds0)
            traj0, t0 = trajectory(ds0, T, u0_T5; Ttr=Ttr, Δt=0.001)
            M0 = Matrix(traj0)
            plot_projection!(pl0,traj0; xlims=(-30, 30), ylims=(-4, 4),
                title=L"Tsit5,dt=%$dt")

            #Rodas5
            u0_R5 = current_state(ds)
            traj, t = trajectory(ds, T, u0_R5; Ttr=Ttr, Δt=0.01)
            plot_projection!(pl,traj; xlims=(-30, 30), ylims=(-4, 4),
                title=L"Rodas5,dt=%$dt")
        end
        push!(plots_adapt, pl, pl0)
    end
    push!(plots,plots_adapt)
end

#save test
#=

begin 
    ds1 = create_ds(u0,p;diffeq=diffeq1)
    ds1_nona = create_ds(u0,p;diffeq=diffeq1_nona)
    ds1_hp = create_ds(u0,p;diffeq=diffeq1_hp)
    ds1_hp_nona = create_ds(u0,p;diffeq=diffeq1_hp_nona)
    ds2 = create_ds(u0,p;diffeq=diffeq2)
    ds2_nona = create_ds(u0,p;diffeq=diffeq2_nona)
    ds2_hp = create_ds(u0,p;diffeq=diffeq2_hp)
    ds2_hp_nona = create_ds(u0,p;diffeq=diffeq2_hp_nona)
end


begin
    traj1_nona,t = trajectory(ds1_nona,T,u0;Ttr=Ttr,Δt=1e-3)
    traj1_hp,t = trajectory(ds1_hp,T,u0;Ttr=Ttr,Δt=1e-3)
    traj1_hp_nona,t = trajectory(ds1_hp_nona,T,u0;Ttr=Ttr,Δt=1e-3)
    traj2_nona,t = trajectory(ds2_nona,T,u0;Ttr=Ttr,Δt=1e-3)
    traj2_hp,t = trajectory(ds2_hp,T,u0;Ttr=Ttr,Δt=1e-3)
    traj2_hp_nona,t = trajectory(ds2_hp_nona,T,u0;Ttr=Ttr,Δt=1e-3)
    traj2,t = trajectory(ds2,T,u0;Ttr=Ttr,Δt=1e-3)
end

begin
    plot_projection!(pl,traj1_nona;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Tsit5 nona, "*L"dt=1e-3")
    plot_projection!(pl,traj1_hp;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Tsit5, "*L"dt=1e-4")
    plot_projection!(pl,traj1_hp_nona;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Tsit5 nona, "*L"dt=1e-4")
    #----
    plot_projection!(pl,traj2_nona;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Rodas5 nona, "*L"dt=1e-3")
    plot_projection!(pl,traj2_hp;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Rodas5, "*L"dt=1e-4")
    plot_projection!(pl,traj2_hp_nona;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Rodas5 nona, "*L"dt=1e-4")
    plot_projection!(pl,traj2;xlims=(-30,30),ylims=(-2,2),la=range(0,1,length(traj1)),label="Rodas5, "*L"dt=1e-3")
end


=#


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
