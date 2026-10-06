using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
cd(@__DIR__)

#--------------IDEA--------------------
#calc lyap curves in function of w, in several parallelized batches 
#take parameter values from fitting of data sets

#system sizes 
N_books = 300
N_mov = 100
N_trends = 100
N_reddit = 100
N_pub = 300
N_wiki = 100
N_twitter = 300 

#constant parameter values for each dataset 
p_const_mov = [25,25,1.0,2.4,0.001,N_mov]
p_const_trends = [22,22,1.0,2.0,0.0005,N_trends]
p_const_reddit = [30,30,1.0, 5.4,0.001,N_reddit]
p_const_pub = [6,6,1.0, 2.4,0.003,N_pub]
p_const_wiki = [48,48,1.0, 2.4,0.0001,N_wiki]

#(w,g) values for each dataset 
w_mov,g_mov = rescaled_params(p_const_mov...)
w_tr,g_tr = rescaled_params(p_const_trends...)
w_red,g_red = rescaled_params(p_const_reddit...)
w_pub,g_pub = rescaled_params(p_const_pub...)
w_wiki,g_wiki = rescaled_params(p_const_wiki...)

#yvalues for each dataset (g_rel_levels)
y_mov = g_mov / (N_mov-1)
y_tr = g_tr / (N_trends-1)
y_red = g_red / (N_reddit-1)
y_pub = g_pub / (N_pub-1)
y_wiki = g_wiki / (N_wiki-1)

#calc time, transient time etc.
T_WLC = 100 #1000
T_fp = 100
Ttr_WLC = 0 #1000
Ttr_fp = 10
#Ttr = 500 #500 #10000 
#T = 500 #500 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9
calc_xsys = false
wstart = 0.01 #[0.65,0.55,0.47] #[0.85,0.8,0.65,0.55,0.47] #0.2
wstop = 0.95 #[0.8,0.75,0.61]  #[0.95,0.92,0.8,0.75,0.61] #0.95
nr_ws = 20 #301 #201
Δw = 0.1
skip =:none #:chaos :fp :none

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)


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

#get lyap curve for a given vertical level on phase diagram
function lyapcurve(p,nr_ws;skip=:none,savedir,T_WLC=1e3,Ttr_WLC=1e4,Ttr_fp=1e2,T_fp=10,log_sys=false,diffeq=diffeq)

    r,_,K,c_ref,α,N = p
    N = Int.(N)

    w_ref,g = rescaled_params(p...)
    g_rel = g/(N-1) #same as y

    #calc critical params
    w_c1 = 1/g_rel #HOPF 
    w_c2 = (N - 1) / N #heteroclinic

    #take w values around w_ref
    w_values = range(w_c1,w_c2;length=nr_ws)

    #set param sweep intervals
    c_start = c_convert.(wstart,r,K,α,N)
    c_stop = c_convert.(w_values[end],r,K,α,N)
    c_values = c_convert.(w_values,r,K,α,N) #same for all system sizes

    mkpath(savedir)
    
    
    Threads.@threads for j in 1:length(c_values) 
        w = w_values[j] 
        c = c_values[j] 
        
        id = Threads.threadid() 
        @show c,id

        #define ds for each thread 
        u0_thread = rand(2N) .* 1e-2 #without the scaling, IT DOESN'T WORK
        p_new = [r, r, K, c, α, N]
        if log_sys
            u0_log = [log.(u0_thread)[1:2:end]...,u0_thread[2:2:end]...] #take log of Ls and reorder
            p_new_log = Union{Int64,Float64}[p[1:end-1]...,Int(N)]
            ds = CoupledODEs(log_attention_in_place!, u0_log, p_new_log;diffeq=diffeq)
        else
            ds = create_original_ds(u0_thread,p_new; diffeq=diffeq) 
        end

        #select calc time and transient time based on w
        Ttr = w_c1 < w < w_c2 ? Ttr_WLC : Ttr_fp
        Ttr_scaled = Ttr / α
        T = w_c1 < w < w_c2 ? T_WLC : T_fp
        T_scaled = T / α
        
        #skip parts of the phase diagram
        if skip==:chaos
            w > w_c1 && continue
        elseif skip==:fp
            w < w_c1 && continue
        end
        
        #calc lyapunov exponent
        #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
        lambda = lyapunov(ds,T_scaled;Ttr=Ttr_scaled)
    
        #lock file access before writing 
        file_lock = ReentrantLock()
    
        lock(file_lock) do
            #open file here
            f = open(savedir*"lyaps_r$(r)_K$(K)_alpha$(α)_N$(N)"*".txt", "a")   
            #save lyapunov exponents
            write(f, "$(g_rel)\t$(r)\t$(w)\t$(c)\t$(lambda)\n")
            close(f) #close file here
        end
    end #of threaded param loop
end


#call for every dataset 
savedir_root = "../../data/lyapcurves_datasets/"
#movies

savedir_mov = savedir_root * "lyapunovs_mov_r$(p_const_mov[1])/"
@time lyapcurve(p_const_mov,nr_ws;
    skip=skip,savedir = savedir_mov,
    T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,
    diffeq=diffeq)
