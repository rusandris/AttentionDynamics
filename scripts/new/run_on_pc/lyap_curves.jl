using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
#using Plots,LaTeXStrings
cd(@__DIR__)

#---------------------------saving-------------------------------

savedir = "./data/data_ensembles/"
mkpath(savedir)

#-----------------------time parameters--------------------------
Tλ_benettin = 2000
Tλ_ens = 10
T_init_traj = 100000
Ttr_init_traj = 10000
Ttrλ_ens = 100.0

@show Tλ_ens 
@show T_init_traj 
@show Ttr_init_traj
#@show Ttrλ_ens 

Ttrsλ_ens = [0,100,500,1000,5000]

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

#test how dependent are the results from T_init_traj

#T_init_traj_values = [100000] #[100,200,500,1000,2000]
λ_avgs_all = Float64[]
for Ttrλ_ens in Ttrsλ_ens
    @show Ttrλ_ens
    λ_avgs = ensemble_average_lyapunov(ds,ws,Tλ_ens,Ttrλ_ens,nr_ens;u0=u0,
        ensemble_init_param=w_ens_init,T_init_traj = T_init_traj,Ttr_init_traj = Ttr_init_traj,Δt = Δt_init_traj_save,filename=savedir*"lyaps_all_Ttrlambda_$Ttrλ_ens.txt")
    
    #@show λ_avgs[1]
    #push!(λ_avgs_all,λ_avgs[1])
    file_name = "lyaps_ensemble_average_T_init_traj$T_init_traj"*"_ens_$nr_ens"*"_Tlambda_$Tλ_ens"* "Ttrlambda_$Ttrλ_ens" *"_N$N"*"_g$g"*"test.txt"
    writedlm(savedir*file_name,hcat(ws,λ_avgs)) 
end


#=
#--------------check if Tλ is enough to calc lyap exps----------------
#random initial conditions
random_u0(dim) = rand(dim) .* 10

Tλs = [100,500,1000,2000,5000]	
for Tλ in Tλs
	@show Tλ 
	λs = calc_lyapunov_curves_reinit(ds,ws,[10000],Tλ;
	generate_u0 = random_u0,d0=1e-9)

	filename_reinit = savedir*"lyapunovs_reinit_g_$g"*"_N_$N"* "Ttr_1e4" * "_T_$Tλ" * ".txt"
	writedlm(filename_reinit,hcat(ws,stack(λs)) ) 
end
=#


#@time calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0,save_dir=savedir)



