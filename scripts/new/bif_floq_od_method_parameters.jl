cd(@__DIR__)
#-----------method parameters------------------
const Ttr = 20000 #transient time 
const T = 5000 #lyapunov time 
const w_min_ratio = 1.01 #w_min/wc_1 for continuation
const w_max_ratio = 0.85 #w_max/wc_2 for continuation


Ns = [3:7;] #nr of competitors dim = 2N 

#select g value from the phase diagram
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
const gval_phase_diag = gs[30]
global_params = [Ns...,gval_phase_diag]

#save params that apply to every system size
#writedlm(cont_dir*"global_params_Ns_g.txt",global_params)

#-----------initial conditions-----------------
#define well-controlled initial condition for each system size

function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end