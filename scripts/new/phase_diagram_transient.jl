using Revise
using DynamicalSystems
using AttentionDynamics
using DelimitedFiles
using OrdinaryDiffEq
cd(@__DIR__)


#=
#read presaved parameters 
data_dir = "../../data/phasediag_data_N10_zoomed/"
saved_files = readdir(data_dir)
params = readdlm(data_dir * "phasediag_params.txt")
ws = params[1:end, 1] 
gs = params[1:end, 2]
=#

data_dir = "../../data/phasediag_data_N10_extended/"
n_param = 60
ws = range(0.5, 0.99; length=n_param)
gs = range(10.8, 27.0; length=n_param)

#diffeq = (alg=Tsit5(), dt=1e-2, adaptive=true, maxiters=1e15)
diffeq2 = (alg=Tsit5(), dt=1e-2, dtmax=1e-3 ,adaptive=true, maxiters=1e15)

#----------------------systems----------------------
N = 10
g = 18.0
w = 0.7
u0 = rand(2N) * 0.1

ds = create_ds(u0, [g, w, N]; diffeq=diffeq2) #no noise

#---------------------time params----------------------
T_fp = 100
Ttr_fp = 100
Ttr_init = 3000
Δt_fp = 0.1
fp_settings = [T_fp, Ttr_fp, Δt_fp]
#Ttr = 2000
Ttrs = [1000] #[100, 500, 1000, 2000] #[100,500,1000,2000]
T = 1000

#lyap settings 
d0 = 1e-9
ut = 1e-6
lt = 1e-12

lyap_vals = Matrix{Float64}[]

#----------------------------------calc lyapunov exponents------------------------------------
@time for (i, Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds, gs, ws; T=T, Ttr=Ttr, Δt=1.0, fp_settings=fp_settings,
        d0=d0, ut=ut, lt=lt, save=false, realtime_plotting=false, parallelize=false)

    push!(lyap_vals, Λ)
end

savedir = "../../data/phasediag_data_N10_extended_d0_$(d0)_ut1e-6/"
mkpath(savedir)
writedlm(savedir * "phasediag_params.txt", hcat(ws, gs))
writedlm(savedir * "lyap_values.txt", lyap_vals[1])
