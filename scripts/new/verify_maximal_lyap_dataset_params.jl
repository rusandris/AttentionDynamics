using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
using CSV,DataFrames
cd(@__DIR__)
#include("load_dataset_param_examples.jl")

#--------------IDEA--------------------
#calc lyap curves in function of w, in several parallelized batches 
#take parameter values from fitting of data sets
dataset_params_dir = "../../data/lyapcurve_datasets_parameters/"
dataset_params_df = CSV.read(dataset_params_dir * "datasets_parameters.txt",DataFrame)

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
w_stop_N100 = 0.1 #Google,Movies,Reddit
w_stop_N300 = Inf #all other datasets
skip =:none #:chaos :fp :none

diffeq = (alg=Tsit5(), dt=dt, dtmax=dtmax, adaptive=true, maxiters=1e15)

#call for every dataset 
savedir_root = "../../data/lyapcurves_datasets_allparams_T$(T_WLC)_res$(nr_ws)/"
mkpath(savedir_root)
#constant params 
K = 1.0
nr_dataset_params = length(dataset_params_df.Dataset)
for i in 1:nr_dataset_params
    dataset_name = dataset_params_df[i,1]
    r,c,α,N,w,y = dataset_params_df[i,2:end]
    params = [r,r,K,c,α,N]         
    w_stop = N == 300 ? w_stop_N300 : w_stop_N100

    @info "Starting $(dataset_name) with $(params), w_stop = $(w_stop)..."
    flush(stderr)

    lyapcurve(params,nr_ws;
        w_stop=w_stop,
        skip=skip,savedir = savedir_root,
        fname_addition="_$(dataset_name)",
        T_WLC=T_WLC,Ttr_WLC=Ttr_WLC,Ttr_fp=Ttr_fp,
        log_sys=true,diffeq=diffeq)
end


