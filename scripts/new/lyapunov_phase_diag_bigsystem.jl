using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots, LaTeXStrings
using LinearAlgebra
cd(@__DIR__)

#---------------------------IDEA-------------------------------
#calc phase diag for N=100 (ORIGINAL EQUATIONS)
#actually change rp and c instead of g and w

data_dir = "../../data/phasediag_data_N10_zoomed/"
#mkpath(data_dir)
figs_dir = "../../paper_figs/"
#mkpath(figs_dir)
saved_files = readdir(data_dir)
params = readdlm(data_dir * "phasediag_params.txt")
ws = params[1:end, 1] #cut some off the end 
gs_N10 = params[1:end, 2]

N = 100
#N = 300
gs = gs_N10 ./ (10-1) .* (N-1)

#original params
#rp=3.0;rc = 3.0;K = 1.0;α = 5*1e-3;c=1.8
r_const = 1.0
rc = r_const
rps = range(6.0,15.0;length=50)
K = 1.0;
cs = range(1.2, 35.0; length=50) 
α = 0.05;
rescaled_params.(rps, rc, K, c, α, N)
rescaled_params.(1, rc, K, cs, α, N)

p = [rps[1], rc, K, cs[1], α]

#critical params from the rescaled 
wc1 = (N - 1) / (rps[1] / α)
wc2 = (N - 1) / N


#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

#create ds

#u01 = 0.1 .* rand(2N)
#u01 = BigFloat.(u0;precision=56) #turn to BigFloat (higher precision)
x0 = rand(N) * 0.001
u02 = zeros(2N)
u02[1:2:end] .= x0

ds = create_original_ds(u02, deepcopy(p); diffeq=diffeq)

T_fp = 100
Ttr_fp = 100
Ttr_init = 3000
Δt_fp = 0.1
fp_settings = [T_fp, Ttr_fp, Δt_fp]
Ttr = 1000
Ttrs = [Ttr]
T = 2000
lyap_vals = zeros(length(rps[1:2:end]),length(cs[1:2:end]))

@time for (i, rp) in enumerate(rps[1:2:end])
    @show i,rp
    wc2 = (N - 1) / N
    ps = current_parameters(ds)
    wc1 = (N - 1) / (rp/α)
    w,g = rescaled_params(rp, rc, K, c, α, N)

    set_parameter!(ds,1,rp)
    for (j,c) in enumerate(cs[1:2:end])
        @show j,c
        set_parameter!(ds,4,c)

        if 0.0 < w < wc1
            λ = lyapunov(ds, T_fp; Ttr=Ttr_fp, d0_lower=1e-12, d0_upper=1e-6, d0=1e-9, Δt=Δt_fp)
        else 
            λ = lyapunov(ds,T;Ttr=Ttr,d0_lower=1e-12,d0_upper=1e-6,d0=1e-9,Δt = 1.0)
        end
        @show λ
        lyap_vals[i,j] = λ
    end
end

mkpath("../../data/big_system_N100_lyapunovs/")
writedlm("../../data/big_system_N100_lyapunovs/lyap_vals.txt",lyap_vals)
writedlm("../../data/big_system_N100_lyapunovs/other_params.txt", (r_const,K,α))
writedlm("../../data/big_system_N100_lyapunovs/rp_c.txt", hcat(rps,cs))

#@time Λ = lyapunov_diagram(ds_noise, gs, ws; T=T, Ttr=Ttr, Δt=1.0, fp_settings=fp_settings,
#    d0=d0, ut=1e-3, lt=1e-14, save=false, realtime_plotting=false, parallelize=false)

#@time Λ = lyapunov_diagram(ds,rps[1:2],ws[1:2]; pidx=[1,4], T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
#d0=1e-9,save=false,realtime_plotting=false,parallelize=false)