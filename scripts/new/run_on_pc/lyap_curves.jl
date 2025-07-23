using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)

#calculate lyapunov curves for `ws` values of w 
#using several different Ttr transient times 
#for the same g 
#without reinit 
function calc_lyapunov_curves(ds,ws,Ttrs,T;u0=current_state(u0))

    lyap_curves = [] #save lyaps here

    for Ttr in Ttrs
        #lyapunov exps of the whole continuous system 
        lyaps = zeros(length(ws))

        #create pds 
        pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9 .* randn(2N)])

        for (i,w) in enumerate(ws)
            @show w
            set_parameter!(pds,2,w)
            lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
        end
        push!(lyap_curves,lyaps)
    end
    return lyap_curves
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
        lyaps_reinit = zeros(length(ws))

        Threads.@threads for i in 1:length(ws)
            w = ws[i]
            @show w
            id = Threads.threadid()
            @show id
            pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9 .* randn(2N)])
            set_parameter!(pds,2,w)
            lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
        end
        push!(lyap_curves,lyaps)
    end
    return lyap_curves
end

#wc1 = (N-1)/g
T = 200
w_min = 0.21
w_max = 0.4
nw = 20
ws = range(w_max,w_min;length=nw) 


diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
p = [g,w_max,N]
u0 = rand(2N) .* 10
#u0 = generate_initial_condition(N)

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)


Ttrs = [10,20,30]
calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0)
calc_lyapunov_curves_reinit(ds,ws,Ttrs,T;u0=u0)