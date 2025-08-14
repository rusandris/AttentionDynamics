using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
#using Plots,LaTeXStrings
cd(@__DIR__)

#---------------------------saving-------------------------------

savedir = "./data/data_benettin/"
mkpath(savedir)

#-----------------------time parameters--------------------------
Tλ_benettin = 2000
Ttrλ_benettin = 100.0
Tλ_ens = 10
T_init_traj = 100000
Ttr_init_traj = 10000
Ttrλ_ens = 100.0

@show Tλ_ens 
@show T_init_traj 
@show Ttr_init_traj
#@show Ttrλ_ens 

Ttrsλ = [0,100,500,1000,5000]

#---------------------ensemble parameters-----------------------
nr_ens = 10000
w_ens_init = 0.3
Δt_init_traj_save = 1.0

@show nr_ens
@show w_ens_init 
@show Δt_init_traj_save


#----------------------system parameters----------------------
#w_min = 0.21
w_min = 0.3
w_max = 0.4
nw = 50 #100
ws = range(w_max,w_min;length=nw) 
#ws = [0.4]
#@show ws
#-----------------------setup for dynamical system----------------------
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10 #N = 5
g = 45.0 #g = 20.0
p = [g,w_max,N]
u0 = generate_initial_condition(N) #dummy init as [1,0,2,0,3...]
u0 = rand(2N) #dummy init close to origin
ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)


#--------------check if TTtr matters when to calc lyap exps----------------
#generate one random initial conditions
#use this as a starting point for every parameter
u0 = rand(2N) 
const_u0(dim) = u0

#Ttrsλ = [0,100,500,1000,5000]
Ttrsλ = [100000]
λs = calc_lyapunov_curves_reinit(ds,ws,Ttrsλ,Tλ_benettin;
	generate_u0 = const_u0,d0=1e-9)

filename_reinit = savedir*"lyapunovs_benettin_long_const_u0_" * "Ttr_lambda_10e5" *"T_$Tλ_benettin" * "_N_$N" * "g_$g" * "_test2.txt"
writedlm(filename_reinit,hcat(ws,stack(λs)) ) 


#@time calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0,save_dir=savedir)



