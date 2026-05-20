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
#Ns = [3,5,8,10,15,20] #system sizes
#Ns = [30,50,75,100]
Ns = [200,300]
g_rel = 2000/299 #4.0 # 1/(\epsilon * (N-1)) is the same for all systems 
gs = g_rel .* (Ns .- 1)
w_c1 = 1/g_rel
ws_c2 = (Ns .- 1) ./ Ns
w_values = [w_c1-0.2:0.002:0.6;]#ws_c2[end]+0.02;]

#original parameter values 
K = 1.0;r=10;α0 = 0.005;
p_L = [r, r, K, 2.4, α0, 300]
c_values = [1.5:0.01:8.5;]

#calc time, transient time etc.
T_WLC = 1000
T_fp = 10
Ttr_WLC = 1000
Ttr_fp = 10
Ttr = 10000
T = 1000
dtmax = 1e-1
dt = 1e-1
diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)

#paths, others
savedir_root = "../../../../../../data/RESULTS/USERS/rusandris/lyapunovs_effect_size/"
savedir = savedir_root * "lyapunovs_parallel_effect_size_Lsys_cs701_largestN_gn6.88_dt1e-1_T$(T)_Ttr$(Ttr)_dtmax_$(dtmax)/"
mkpath(savedir)
writedlm(savedir * "system_sizes.txt", Ns)

Threads.@threads for i in 1:length(Ns)
    N = Ns[i]
    id = Threads.threadid()
    @info "N = $(N), thread id = $(id)"
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    #g = gs[i]
    #wc2 = ws_c2[i]

    #generate standard init cond
    u0 = rand(2N)
    #u0 = generate_initial_condition(N)
    #@show u0

    #define system 
    #p = [g,w_values[end],N]
    #create system
    #ds = create_ds(u0,p;diffeq=diffeq)
    p = [p_L...,N]
    ds = create_original_ds(u0, p; diffeq=diffeq)

    #create pds 
    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9])

    #array for lyapunov exps 
    lyaps = zeros(length(w_values))
    
    for (j,c) in enumerate(c_values)
        #@show w
        @info "w = $(c), thread id = $(id)"
	reinit!(pds)
        #set_parameter!(pds,2,w)

	set_parameter!(pds,4,c)

        #select calc time and transient time based on w
        #Ttr = w_c1 < w < wc2 ? Ttr_WLC : Ttr_fp
        #T = w_c1 < w < wc2 ? T_WLC : T_fp

        #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
	lambda = lyapunov(pds,T;Ttr=Ttr)
	open(savedir*"lyaps_N$N"*".txt", "a") do f
		write(f, "$(c)\t$(lambda)\n")
	end

    end  

    #writedlm(savedir*"lyaps_N$N"*".txt",hcat(w_values,lyaps))
end
