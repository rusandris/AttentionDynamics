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
T_WLC = 500 #1000
T_fp = 500
Ttr_WLC = 0 #1000
Ttr_fp = 10
#Ttr = 500 #500 #10000 
#T = 500 #500 #1000
dtmax = 1e-2
dt = 1e-2
d0 = 1e-9
calc_xsys = false
nr_ws = 100 #301 #201
Δw = 0.1
skip =:none #:chaos :fp :none

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)


#call for every dataset 
savedir_root = "../../data/lyapcurves_datasets/"
#movies

savedir_mov = savedir_root * "lyapunovs_mov_test_r$(p_const_mov[1])/"
@time lyapcurve(p_const_mov,nr_ws;
    skip=skip,savedir = savedir_mov,
    T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,
    diffeq=diffeq)

savedir_tr = savedir_root * "lyapunovs_tr_r$(p_const_trends[1])/"
@time lyapcurve(p_const_trends,nr_ws;
    skip=skip,savedir = savedir_tr,
    T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,
    diffeq=diffeq)

savedir_red = savedir_root * "lyapunovs_red_r$(p_const_reddit[1])/"
@time lyapcurve(p_const_reddit,nr_ws;
    skip=skip,savedir = savedir_red,
    T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,
    diffeq=diffeq)

savedir_pub = savedir_root * "lyapunovs_pub_r$(p_const_pub[1])/"
@time lyapcurve(p_const_pub,nr_ws;
    skip=skip,savedir = savedir_pub,
    T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,log_sys=true,
    diffeq=diffeq)
