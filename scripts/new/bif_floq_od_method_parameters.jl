cd(@__DIR__)
#-----------method parameters------------------
Ttr = 20000 #transient time 
T = 5000 #lyapunov time 
w_min_ratio = 0.9 #w_min/wc_1 for continuation
w_max_ratio = 0.6 #w_max/wc_2 for continuation
w_min_ratio_od = 1.005 #w_min/wc_1 for orbit diag

#Ns = [2:6;] #nr of competitors dim = 2N 
#Ns = [3,5,10] #[2] #[2,5,10]#[10] #
Ns = [2] 

#select g value from the phase diagram
#gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
#gs = gs_N10 ./(10-1)
#const gval_phase_diag = gs[30]
#gval_phase_diag = 5.0
g_rel = 5.0
global_params = [Ns...,g_rel]


#save params that apply to every system size
#writedlm(cont_dir*"global_params_Ns_g.txt",global_params)

