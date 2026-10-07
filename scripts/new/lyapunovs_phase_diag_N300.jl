using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
cd(@__DIR__)

#--------------IDEA--------------------
#calc entire phase diagram for N=300, in several parallelized batches 

Ns = [300] #[20,30,50,75] #[10,200,300] #[100,200,300]
N = Ns[1]
#read g_rel values for this batch from file given as CL ARG
batch_source = ARGS[1]
g_rel_values = readdlm(batch_source)

#read B2,B4 levels from pre-generated file
#params_phase_diag_wg = readdlm("../../data/parameters/alpha_rp=rc_article/params_w_g.txt")
#g_rel_values = [params_phase_diag_wg[:,end][end]] 

#original parameter values 
K = 1.0; α0 = 0.005;
p_L_const = [K,α0] #constant parameters (L system)
#calc time, transient time etc.
T_WLC = 500 #1000
T_fp = 10
Ttr_WLC = 500 #1000
Ttr_fp = 10
#Ttr = 500 #500 #10000 
#T = 500 #500 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9
calc_xsys = false
wstart = 0.01 #[0.65,0.55,0.47] #[0.85,0.8,0.65,0.55,0.47] #0.2
wstop = 0.95 #[0.8,0.75,0.61]  #[0.95,0.92,0.8,0.75,0.61] #0.95
nr_ws = 80 #301 #201
skip =:none #:chaos :fp :none

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
diffeqx = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)

#diffeq = (alg=Rodas5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)
#diffeqx = (alg=Rodas5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)
#paths, others
savedir_root = "../../data/phase_diag_N300_recalc_fixed/" #"../../../../../../data/RESULTS/USERS/rusandris/lyapunovs_effect_size/"
 
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


#function call
@time phasediag_lyapcurves(Ns,g_rel_values;wstart=wstart,wstop=wstop,nr_ws=nr_ws,skip=skip,T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,p_L_const = p_L_const,diffeq=diffeq,diffeqx=diffeqx)

