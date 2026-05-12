using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
cd(@__DIR__)

#--------------IDEA--------------------
#show that width of chaotic band grows with system size N 

#system sizes, parameters
Ns = [3,5,8,10,15,20] #system sizes 
g_rel = 4.0 # 1/(\epsilon * (N-1)) is the same for all systems 
gs = g_rel .* (Ns .- 1)
w_c1 = 1/g_rel
ws_c2 = (Ns .- 1) ./ Ns
w_values = [w_c1:0.1:ws_c2[end];]

#calc time, transient time etc.
T = 1000
Ttr_WLC = 10000
Ttr_fp = 10 
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e15)

#paths, others 
savedir = "../../data/lyapunovs_effect_size/"
mkpath(savedir)

for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = gs[i]
    wc2 = ws_c2[i]

    #generate standard init cond
    u0 = rand(2N)
    #u0 = generate_initial_condition(N)
    #@show u0

    #define system 
    p = [g,w_values[end],N]
    #create system
    ds = create_ds(u0,p;diffeq=diffeq)

    #create pds 
    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9])

    #array for lyapunov exps 
    lyaps = zeros(length(w_values))

    for (i,w) in enumerate(w_values)
        @show w
        set_parameter!(pds,2,w)

        #select transient time based on w
        Ttr =  w < wc2 ? Ttr_WLC : Ttr_fp
        lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
    end  

    writedlm(savedir*"lyaps_N$N"*".txt",hcat(w_values,lyaps))
end