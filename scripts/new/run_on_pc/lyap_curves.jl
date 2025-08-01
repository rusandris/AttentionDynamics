using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
#using Plots,LaTeXStrings
cd(@__DIR__)

#---------------------------saving-------------------------------

savedir = "./data_highres/"
mkpath(savedir)

#-----------------------time parameters--------------------------
T = 200
Tλ = 50
T_init_traj = 2000
Ttr_init_traj = 5000
Ttrs = [100,500,1000,5000,10000,500000,100000]


#----------------------system parameters----------------------
w_min = 0.21
w_max = 0.4
nw = 10
ws = range(w_max,w_min;length=nw) 

#-----------------------setup for dynamical system----------------------
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10 #N = 5
g = 45.0 #g = 20.0
p = [g,w_max,N]
u0 = generate_initial_condition(N) #dummy init as [1,0,2,0,3...]
ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)

#random initial conditions
random_u0(dim) = rand(dim) .* 10

λs =  calc_lyapunov_curves_reinit(ds,ws,[1000],2000;
        generate_u0 = random_u0,d0=1e-9)

filename_reinit = savedir*"lyapunovs_reinit_g_$g"*"_N_$N"*".txt"
writedlm(filename_reinit,hcat(ws,stack(λs)) ) 

#@time calc_lyapunov_curves(ds,ws,Ttrs,T;u0=u0,save_dir=savedir)


#test ensemble lyap calc
#=
ens_res = [] 
for Ttr in Ttrs[1:4]
    @show Ttr
    λ_ens = ensemble_average_lyapunov(ds,ws,T,Ttr,10;u0=current_state(ds),
        ensemble_init_param=0.3,T_init_traj = T,Ttr_init_traj = 1000,Δt = 0.01)
    push!(ens_res,λ_ens)
end
=#