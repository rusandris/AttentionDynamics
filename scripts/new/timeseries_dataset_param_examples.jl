using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using OrdinaryDiffEq
cd(@__DIR__)

include("load_dataset_param_examples.jl")

labels = ["Movies 1995-2000","Publications 1990-2010","Reddit 2012","Google Trends 2014","Books 1930-1950"]
diffeq = (alg=Tsit5(), dt=1e-2, dtmax=1e-2, adaptive=true, maxiters=1e15)

#times in x system
Tx = 500
Ttrx = 500

#movies 
u0 = rand(2N_mov) .* 1e-2
u0_mov = [log.(u0)[1:2:end]...,u0[2:2:end]...]
p_log_mov = Union{Int64,Float64}[p_const_mov[1:end-1]...,Int(p_const_mov[end])]
ds_mov = create_logtransform_ds(u0_mov,p_log_mov;diffeq=diffeq) #mov

α_mov = p_const_mov[5]

T_mov = 100 #Tx / α_mov
Ttr_mov = 5000 #Ttrx / α_mov

#dataset param 
traj_mov,t = trajectory(ds_mov, T_mov; Ttr = Ttr_mov)
lt = Matrix(traj_mov)[:,1:N_mov]
Lt = exp.(lt)

plot(t,Lt,legend=false,xlabel=L"t",ylabel=L"L(t)",
    guidefontsize=18,tickfontsize=18,title = L"y=%$(y_mov),w=%$(w_mov)")

#move to the right (increase w)
p_mod_mov = Union{Int64,Float64}[p_const_mov[1:3]...,84.0,α_mov,Int(N_mov)]
rescaled_params(p_mod_mov...) 
ds_mov = create_logtransform_ds(u0_mov,p_mod_mov;diffeq=diffeq) #mov
traj_mov_mod,t = trajectory(ds_mov, T_mov; Ttr = Ttr_mov)
lt = Matrix(traj_mov_mod)[:,1:N_mov]
Lt = exp.(lt)

plot(t,Lt,legend=false,xlabel=L"t",ylabel=L"L(t)",
    guidefontsize=18,tickfontsize=18,title = L"y=%$(y_mov),w=%$(w_mov)")


#plot and save each selected example param from the fitted datasets 
dataset_params = [p_const_mov,p_const_pub,p_const_reddit,p_const_trends,p_const_books]
T = 100 #in LY time 
Ttr = 5000 ##in LY time
plots_ts = []
for (i,params) in enumerate(dataset_params)

    @show labels[i]
    @show params
    α = params[5]
    N = Int(params[end])

    u0 = rand(2N) .* 1e-2
    u0_log = [log.(u0)[1:2:end]...,u0[2:2:end]...]
    p_log = Union{Int64,Float64}[params[1:end-1]...,Int(params[end])]
    ds_log = create_logtransform_ds(u0_log,p_log;diffeq=diffeq) 

 

    #dataset param 
    traj_log,t = trajectory(ds_log, T; Ttr = Ttr)
    lt = Matrix(traj_log)[:,1:N]
    Lt = exp.(lt)

    pl = plot(t,Lt,legend=false,xlabel=L"t",ylabel=L"L(t)",title = labels[i])

    push!(plots_ts,pl)
end

pl_ts = plot(plots_ts...,layout = (length(plots_ts),1),
    size=(800,1800),guidefontsize=18,tickfontsize=10,
    left_margin=15Plots.mm)

savefig(pl_ts,"../../paper_figs/dataset_param_examples_timeseries_LY_logsys_trans_T$(T)_Ttr_$(Ttr).png")