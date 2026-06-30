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
Ns = [300] #[20,30,50,75] #[10,200,300] #[100,200,300]
#g_rel_values = [4.456843049327354,2000/299] #[2.0, 3.0, 4.456843049327354, 5.0, 2000/299] #[2.0,3.0] #4.0 # 1/(\epsilon * (N-1)) is the same for all systems 
#read B2,B4 levels from pre-generated file
params_phase_diag_wg = readdlm("../../data/parameters/alpha_rp=rc_article/params_w_g.txt")
g_rel_values = [params_phase_diag_wg[:,end][end]] 

#original parameter values 
K = 1.0; α0 = 0.005;c0 = 2.4
nr_r = 50
r_values_ends = readdlm("../../data/parameters/alpha_rp=rc_article/params_r.txt")[[1,end]]
r_values = range(r_values_ends[1],r_values_ends[end],length=nr_r)
p_L_const = [K,α0] #constant parameters (L system)
c_values = [2.4]
p_log = Union{Int64,Float64}[10.0,10.0,1.0,2.4,0.005,300]

#calc time, transient time etc.
T_WLC = 1000
T_fp = 10
Ttr_WLC = 1000
Ttr_fp = 10
Ttr = 500 #500 #10000
Ttr_lyap = 10 
T = 500 #500 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9
calc_xsys = false
wstart = 0.15 #[0.65,0.55,0.47] #[0.85,0.8,0.65,0.55,0.47] #0.2
wstop = 0.95 #[0.8,0.75,0.61]  #[0.95,0.92,0.8,0.75,0.61] #0.95
nr_ws = 200 #301 #201

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
diffeqx = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)

#diffeq = (alg=Rodas5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
#diffeqx = (alg=Rodas5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
#paths, others
savedir_root = "../../data/lyapunovs_effect_size_L_N300_logsys_change_r/" #"../../../../../../data/RESULTS/USERS/rusandris/lyapunovs_effect_size/"
#writedlm(savedir * "gs_alphas.txt", hcat(g_rel_values,αs))

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
function lyapunov_exp_curves(Ns,c_values,r_values;T=1e3,Ttr=1e4,Ttr_lyap=1e2,p_L_const = [10,10,0.005],log_sys=false,diffeq=diffeq,diffeqx=diffeqx,calc_xsys=false)
    c_values_shortened = round.(c_values;digits=2) 
    K,α0 = p_L_const
    for (ci,c) in enumerate(c_values)
        c_short = c_values_shortened[ci] 
       
        #scaled calc times, transient times
        Ts_scaled = T ./ α0 
        Ttrs_scaled = Ttr ./ α0
       
        savedir = savedir_root * "lyapunovs_parallel_effect_size_Lsys__$(c_short)_scaled_times_T$(T)_Ttr$(Ttr)_dtmax_$(dtmax)/"
        mkpath(savedir)
        writedlm(savedir * "system_sizes.txt", Ns)
        #writedlm(savedir * "gs_alphas.txt", hcat(gs,αs))
        writedlm(savedir * "times.txt", hcat(Ts_scaled,Ttrs_scaled))

        for i in 1:length(Ns)
            N = Ns[i]
            @info "N = $(N)"

            #select g = 1/ε for each N (system size)
            #so g/(N-1) is the same for each N
            #g = gs[i]
            #generate standard init cond
            u0 = rand(2N)
            #u0 = generate_initial_condition(N)
            #@show u0

            #define system 
            #p = [g,w_values[end],N]
            #ds = create_ds(u0,p;diffeq=diffeq)
            #p = [p_L[1:end-2]...,α,N]
           
           
            Threads.@threads for j in 1:length(r_values)  
                r = r_values[j]
                g_rel_calc = r/(α0*(N-1))
                T_scaled = Ts_scaled[i] 
                Ttr_scaled = Ttrs_scaled[i] 
                id = Threads.threadid() 
                @info "1/eps/(N-1) = $(g_rel_calc),r = $(r),c = $(c), thread id = $(id)" 
                flush(stderr) 
         
                #define ds for each thread 
                u0_thread = rand(2N) .* 1e-2 #without the scaling, IT DOESN'T WORK
                p_new = [r,r,K,c,α0, N]
                if log_sys
                    u0_log = [log.(u0_thread)[1:2:end]...,u0_thread[2:2:end]...] #take log of Ls and reorder
                    p_new_log = Union{Int64,Float64}[p_new[1:end-1]...,Int(N)]
                    ds = CoupledODEs(log_attention_in_place!, u0_log, p_new_log;diffeq=diffeq)
                else
                    ds = create_original_ds(u0_thread,p_new; diffeq=diffeq) 
                end
                #set_parameter!(ds, 5, c) 

                #transient time separately 
                u0_Ttr = trajectory(ds, 0; Ttr=Ttr_scaled)[1][1] 
                null_idx = findfirst(x->iszero(x),u0_Ttr) 
                #@show u0_Ttr[1] 
                #@show null_idx, u0_Ttr[null_idx] 

                flush(stdout)
                flush(stderr)

                #detect trapped trajectories
                istrapped = any(iszero.(u0_Ttr))
                if istrapped
                        @warn "From thread ID $(id): Zeros in the system's state at r = $r, skipping value... "
                    
                #write out blank instead of lyap
                open(savedir*"lyaps_N$N"*".txt", "a") do f
                    write(f, "$(g_rel_calc)\t$(r)\t\n")
                end            
                    continue 
                end

                #create pds 
                pds = ParallelDynamicalSystem(deepcopy(ds),[u0_Ttr,u0_Ttr .+ d0])

                #select calc time and transient time based on w
                #Ttr = w_c1 < w < wc2 ? Ttr_WLC : Ttr_fp
                #T = w_c1 < w < wc2 ? T_WLC : T_fp

                #calc lyapunov exponent
                #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
                lambda = lyapunov(pds,T_scaled;Ttr=Ttr_lyap / α0)
        
                #save lyapunov exponents
                open(savedir*"lyaps_N$N"*".txt", "a") do f
                     write(f, "$(g_rel_calc)\t$(r)\t$(lambda)\n")
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
@time lyapunov_exp_curves(Ns,c_values,r_values;T=T,Ttr=Ttr,Ttr_lyap=Ttr_lyap,log_sys=true,p_L_const = p_L_const,diffeq=diffeq,diffeqx=diffeqx)

