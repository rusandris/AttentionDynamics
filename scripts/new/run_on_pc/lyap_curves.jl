using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
#using Plots,LaTeXStrings
cd(@__DIR__)

#calculate lyapunov curves for `ws` values of w 
#using several different Ttr transient times 
#for the same g 
#without reinit 
function calc_lyapunov_curves(ds,ws,Ttrs,T;u0=current_state(u0),save_dir="./")

    lyap_curves = [] #save lyaps here

    for Ttr in Ttrs
        g,_,N = current_parameters(ds)
        N = Int(N)
        #create pds 
        pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9 .* randn(2N)])
        file_name = "lyapunovs_g_$g"*"_N_$N"*"_Ttr_$Ttr"*".txt"

        for (i,w) in enumerate(ws)
            println("Ttr: $Ttr" * " w = $w")
            set_parameter!(pds,2,w)
            λ = lyapunov(pds,T;Ttr=Ttr)
            
            open(save_dir*file_name, "a") do f
                write(f, "$w \t $λ \n")
            end

        end
    end
end

#calculate lyapunov curves for `ws` values of w 
#using several different Ttr transient times 
#for the same g 
#with reinit
#paralellized version 
function calc_lyapunov_curves_reinit(ds,ws,Ttrs,T;u0=current_state(u0))

    lyap_curves = [] #save lyaps here

    for Ttr in Ttrs
        #lyapunov exps of the whole continuous system 
        lyaps = zeros(length(ws))

        Threads.@threads for i in 1:length(ws)
            w = ws[i]
            id = Threads.threadid()
            println("Id: $id" * " Ttr: $Ttr" * " w = $w")
            pds = ParallelDynamicalSystem(deepcopy(ds),[u0,u0 .+ 1e-9 .* randn(2N)])
            set_parameter!(pds,2,w)
            lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
        end
        push!(lyap_curves,lyaps)
    end
    return lyap_curves
end

#initialize ensemble from the chaotic attractor?
#for each parameter, take the average LLE  
function ensemble_average_lyapunov(ds,ws,T,Ttr,nr_ens;u0=current_state(ds),
        ensemble_init_param,T_init_traj = T,Ttr_init_traj = Ttr,Δt = 0.01)

    lyap_avgs = zeros(length(ws)) #save lyap avgs here
    

    #init ensemble from a trajectory
    set_parameter!(ds,2,ensemble_init_param)
    traj,t = trajectory(ds,T_init_traj;Ttr=Ttr_init_traj)

    for i in 1:length(ws)
        w = ws[i]

        #go over initial conditions (paralellized)
        #select nr_ens init conds from traj
        lyaps_ensemble = zeros(length(ws)) #save lyap avgs here

        Threads.@threads for e in 1:nr_ens
            id = Threads.threadid()
            println("Id: $id" * " Ttr: $Ttr" * " w = $w")

            #initialize from previously computed traj
            ic_idx = rand(1:length(traj))
            u0 = traj[ic_idx]
            pds = ParallelDynamicalSystem(deepcopy(ds),[u0,u0 .+ 1e-9 .* randn(2N)])
            set_parameter!(pds,2,w)
            lyaps_ensemble[e] = lyapunov(pds,T)
        end
        lyap_avgs[i] = mean(lyaps_ensemble)
    end

    return lyap_avgs
end

#wc1 = (N-1)/g
T = 200
w_min = 0.21
w_max = 0.4

nw = 100

ws = range(w_max,w_min;length=nw) 

savedir = "./data_highres/"
mkpath(savedir)

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0

#N = 5
#g = 20.0
p = [g,w_max,N]
#u0 = rand(2N) .* 10
u0 = generate_initial_condition(N)

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)


Ttrs = [100,500,1000,5000,10000,500000,100000]

#@time calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0,save_dir=savedir)
#@time λs = calc_lyapunov_curves_reinit(ds,ws,Ttrs,T;u0=u0)
#filename_reinit = savedir*"lyapunovs_reinit_g_$g"*"_N_$N"*".txt"
#writedlm(filename_reinit,hcat(ws,stack(λs)) ) 

#test ensemble lyap calc 
Ttr = 0.0
λ_ens = ensemble_average_lyapunov(ds,ws,T,Ttr,100;u0=current_state(ds),
        ensemble_init_param=0.3,T_init_traj = T,Ttr_init_traj = 1000,Δt = 0.01)