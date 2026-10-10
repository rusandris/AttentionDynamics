using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
cd(@__DIR__)
include("load_dataset_param_examples.jl")

#--------------IDEA--------------------
#calc lyap curves in function of w, in several parallelized batches 
#take parameter values from fitting of data sets



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
