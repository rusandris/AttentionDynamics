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
Ns = [20,30,50,75] #[100,150,200] #[20,30,50,75] #[50,75] #[300] #[20,30,50,75] #[10,200,300] #[100,200,300]
g_rel_values = [1.5, 2.0, 2.981165919282511, 4.0, 5.0, 6.0, 2000/299, 8.0, 9.0, 10.0] # 1/(\epsilon * (N-1)) is the same for all systems 
#original parameter values 
K = 1.0;r=10;α0 = 0.005;
p_L_const = [r, r, K] #constant parameters (L system)

#calc time, transient time etc.
T_WLC = 1000
T_fp = 10
Ttr_WLC = 1000
Ttr_fp = 10
Ttr = 50000 #50000 #1000 #10000
Ttr_lyap = 100 
T = 1000 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9
calc_xsys = false
nr_ws = 400 #40 #351 #201 #101
wstart = 0.2 #[0.85,0.8,0.6,0.55,0.47] #0.2
wstop = 0.95 #[0.95,0.92,0.8,0.75,0.61] #0.95

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
diffeqx = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)

#diffeq = (alg=Rodas5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
#diffeqx = (alg=Rodas5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
#paths, others
savedir_root = "../../data/lyapunovs_effect_size_L_logsys_fixed_Ttr_5x1e4_0.2_0.95_w400/" #"../../../../../../data/RESULTS/USERS/rusandris/lyapunovs_effect_size/"
 
c_convert(w,rc,K,α,N) = (rc/(α*K*(N-1))) / (1/w - 1) #convert w values to corresponding c values

function log_attention_in_place!(du, u, p, t)
    rp, rc, K, c, α, N = p
    l = @view u[1:N]
    Y = @view u[N+1:2*N]

    du_l = @view du[1:N]
    du_Y = @view du[N+1:2*N]

    total_attention = 0.0
    @inbounds for j in 1:N
        total_attention += exp(l[j])
    end

    @inbounds for i in 1:N
        sum_competitors = total_attention - exp(l[i])
        du_l[i] = rp * (1.0 - (rc / K) * Y[i] - c * sum_competitors)
        du_Y[i] = exp(l[i]) - α * Y[i]
    end
    return nothing
end

#loop over eps params (vertical levels on phase diagram)
function lyapunov_exp_curves(Ns,g_rel_values;wstart=0.5,wstop=0.8,nr_ws=201,T=1e3,Ttr=1e4,Ttr_lyap=1e2,log_sys=false,p_L_const = [10,10,0.005],diffeq=diffeq,diffeqx=diffeqx,calc_xsys=false)
    g_rel_level_shortened = round.(g_rel_values;digits=2) 

    for (gi,g_rel) in enumerate(g_rel_values)
        glevel_short = g_rel_level_shortened[gi] 
        
        #calc params
        gs = g_rel .* (Ns .- 1) 
        w_c1 = 1/g_rel #HOPF 
        ws_c2 = (Ns .- 1) ./ Ns #heteroclinic
        αs = r ./ (g_rel .* (Ns .- 1) ) #compute α for each system size
     
        #scaled calc times, transient times
        Ts_scaled = T ./ αs 
        Ttrs_scaled = Ttr ./ αs
       
        #set param sweep intervals
        wstart_g = length(wstart) == 1 ? wstart : wstart[gi]
        wstop_g = length(wstop) == 1 ? wstop : wstop[gi]
        w_values = range(wstart_g,wstop_g;length=nr_ws) #[w_c1:0.002:0.8;]
        c_start = c_convert.(w_c1,r,K,αs[1],Ns[1])
        c_stop = c_convert.(w_values[end],r,K,αs,Ns)
        #c_values = [1.5:0.01:8.5;]
        c_values = c_convert.(w_values,r,K,αs[1],Ns[1]) #same for all system sizes

        savedir = savedir_root * "lyapunovs_parallel_effect_size_Lsys__$(glevel_short)_scaled_times_T$(T)_Ttr$(Ttr)_dtmax_$(dtmax)/"
        mkpath(savedir)
        writedlm(savedir * "system_sizes.txt", Ns)
        writedlm(savedir * "gs_alphas.txt", hcat(gs,αs))
        writedlm(savedir * "times.txt", hcat(Ts_scaled,Ttrs_scaled))

        for i in 1:length(Ns)
            N = Ns[i]
            @info "N = $(N)"

            #select α
            α = αs[i]

            #generate standard init cond
            u0 = rand(2N)
            #u0 = generate_initial_condition(N)
            #@show u0

            g_rel_calc = r/(α*(N-1))

            #array for lyapunov exps 
            lyaps = zeros(length(w_values))
            
            Threads.@threads for j in 1:length(c_values) 
                w = w_values[j] 
                c = c_values[j] 
                T_scaled = Ts_scaled[i] 
                Ttr_scaled = Ttrs_scaled[i] 
                id = Threads.threadid() 
                @info "1/eps/(N-1) = $(g_rel_calc),w = $(w),c = $(c), thread id = $(id)" 
                flush(stderr) 
      
                #define ds for each thread 
                u0_thread = rand(2N) .* 1e-2 #without the scaling, IT DOESN'T WORK
                p_new = [p_L_const..., c, α, N]
                if log_sys
                    u0_log = [log.(u0_thread)[1:2:end]...,u0_thread[2:2:end]...] #take log of Ls and reorder
                    p_new_log = Union{Int64,Float64}[p_new[1:end-1]...,Int(N)]
                    ds = CoupledODEs(log_attention_in_place!, u0_log, p_new_log;diffeq=diffeq)
                else
                    ds = create_original_ds(u0_thread,p_new; diffeq=diffeq)
                end
                
                flush(stdout)
                flush(stderr)

                #create pds 
                pds = ParallelDynamicalSystem(deepcopy(ds),[u0_thread,u0_thread .+ d0])

                #select calc time and transient time based on w
                #Ttr = w_c1 < w < wc2 ? Ttr_WLC : Ttr_fp
                #T = w_c1 < w < wc2 ? T_WLC : T_fp

                #calc lyapunov exponent
                #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
                lambda = lyapunov(pds,T_scaled;Ttr=Ttr_scaled)
        
                #save lyapunov exponents
                open(savedir*"lyaps_N$N"*".txt", "a") do f
                     write(f, "$(g_rel_calc)\t$(w)\t$(c)\t$(lambda)\n")
                end
                
                if calc_xsys
                    u0_Ttrx = trajectory(dsx, 0; Ttr=Ttr)[1][1] 
                    dsx = create_ds(rand(2N), [g_rel*(N-1),w, N]; diffeq=diffeqx)
                    pdsx = ParallelDynamicalSystem(deepcopy(dsx),[u0_Ttrx,u0_Ttrx .+ d0])
                    lambdax = lyapunov(pdsx,T;Ttr=Ttr_lyap)

                     #save lyapunov exponents
                    open(savedir*"lyapsx_N$N"*".txt", "a") do f
                         write(f, "$(g_rel_calc)\t$(w)\t$(lambdax)\n")
                    end
                end #of x sys calc

            end #of threaded param loop
        end  #of system size
    end #of eps level param
end

#function call
@time lyapunov_exp_curves(Ns,g_rel_values;nr_ws=nr_ws,wstart=wstart,wstop=wstop,log_sys=true,T=T,Ttr=Ttr,Ttr_lyap=Ttr_lyap,p_L_const = p_L_const,diffeq=diffeq,diffeqx=diffeqx)

