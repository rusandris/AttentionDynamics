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

#wc1 = (N-1)/g
T = 2000
w_min = 0.21
w_max = 0.4
<<<<<<< HEAD
nw = 100
=======
nw = 10
>>>>>>> a1b00a4a648ed138c646ee9b0d6e65c73cb53a87
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


<<<<<<< HEAD
Ttrs = [100,500,1000,5000,10000,500000,100000]
=======
Ttrs = [500,1000,2000]
>>>>>>> a1b00a4a648ed138c646ee9b0d6e65c73cb53a87
#@time calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0,save_dir=savedir)
@time λs = calc_lyapunov_curves_reinit(ds,ws,Ttrs,T;u0=u0)
filename_reinit = savedir*"lyapunovs_reinit_g_$g"*"_N_$N"*".txt"
writedlm(filename_reinit,hcat(ws,stack(λs)) ) 
