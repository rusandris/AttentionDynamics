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
Ns = [100,200,300]
g_rel = 2000/299 #4.0 # 1/(\epsilon * (N-1)) is the same for all systems 
gs = g_rel .* (Ns .- 1)
w_c1 = 1/g_rel #HOPF 
ws_c2 = (Ns .- 1) ./ Ns #heteroclinic
w_values = [w_c1:0.002:0.8;]#ws_c2[end]+0.02;]

#original parameter values 
K = 1.0;r=10;α0 = 0.005;
p_L = [r, r, K, 2.4, α0, 300]
αs = r ./ (g_rel .* (Ns .- 1) ) #compute α for each system size

#convert w values to corresponding c values
c_convert(w,rc,K,α,N) = (rc/(α*K*(N-1))) / (1/w - 1)
c_start = c_convert.(w_c1,r,K,αs[1],Ns[1])
c_stop = c_convert.(w_values[end],r,K,αs,Ns)
#c_values = [1.5:0.01:8.5;]
c_values = c_convert.(w_values,r,K,αs[1],Ns[1]) #same for all system sizes

#calc time, transient time etc.
T_WLC = 1000
T_fp = 10
Ttr_WLC = 1000
Ttr_fp = 10
Ttr = 1000 #1000 #10000
T = 500 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9

#scaled calc times, transient times
Ts_scaled = T ./ αs 
Ttrs_scaled = Ttr ./ αs

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)

#paths, others
savedir_root = "../../data/lyapunovs_effect_size/" #"../../../../../../data/RESULTS/USERS/rusandris/lyapunovs_effect_size/"
savedir = savedir_root * "lyapunovs_parallel_effect_size_Lsys_test_cs_largestN_gn6.88_dt1e-1_T$(T)_Ttr$(Ttr)_dtmax_$(dtmax)/"
mkpath(savedir)
writedlm(savedir * "system_sizes.txt", Ns)
writedlm(savedir * "gs_alphas.txt", vcat(gs,αs))

for i in 1:length(Ns)
    N = Ns[i]
    @info "N = $(N)"

    #select g = 1/ε for each N (system size)
    #so g/(N-1) is the same for each N
    #g = gs[i]

    #select α
    α = αs[i]

    #generate standard init cond
    u0 = rand(2N)
    #u0 = generate_initial_condition(N)
    #@show u0

    #define system 
    #p = [g,w_values[end],N]
    #ds = create_ds(u0,p;diffeq=diffeq)
    p = [p_L[1:end-2]...,α,N]

    #array for lyapunov exps 
    lyaps = zeros(length(w_values))
    
    Threads.@threads for j in 1:length(c_values)
        w = w_values[j]
        c = c_values[j]
        id = Threads.threadid()
        @info "w = $(c), thread id = $(id)"
        flush(stderr)

        #define ds for each thread
        u0_thread = rand(2N) .* 1e-2
        ds = create_original_ds(u0_thread, [p_L[1:end-3]..., c, α, N]; diffeq=diffeq)
        #set_parameter!(ds, 4, c)
        
        #transient time separately
        u0_Ttr = trajectory(ds, 0, u0_thread; Ttr=Ttr)[1][1]
        null_idx = findfirst(x->iszero(x),u0_Ttr)
        @show u0_Ttr[1]
        #@show null_idx, u0_Ttr[null_idx]
        flush(stdout)
        flush(stderr)

        #detect trapped trajectories
        istrapped = any(iszero.(u0_Ttr))
        if istrapped
            @warn "From thread ID $(id): Zeros in the system's state at c = $c, aborting... "
            #write out NaN
            open(savedir*"lyaps_N$N"*".txt", "a") do f
		        write(f, "$(c)\tNaN\n")
	        end            
            continue 
        end

        #create pds 
        pds = ParallelDynamicalSystem(deepcopy(ds),[u0_Ttr,u0_Ttr .+ d0])
        #set_parameter!(pds,2,w)

        #select calc time and transient time based on w
        #Ttr = w_c1 < w < wc2 ? Ttr_WLC : Ttr_fp
        #T = w_c1 < w < wc2 ? T_WLC : T_fp

        #calc lyapunov exponent
        #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
	    lambda = lyapunov(pds,T;Ttr=0)

        #save lyapunov exponents
	    open(savedir*"lyaps_N$N"*".txt", "a") do f
		    write(f, "$(w)\t$(c)\t$(lambda)\n")
	    end

    end  

    #writedlm(savedir*"lyaps_N$N"*".txt",hcat(w_values,lyaps))
end
