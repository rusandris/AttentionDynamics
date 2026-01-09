using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")


#N=7 chaotic regime (g,w,N =  [10.051282051282051, 0.65, 7.0]) dissapears for Ttr = 10^7
#N=6 it dissapears for shorter transients

#read w and g params
#ws_saved = readdlm("../../plotting/data/N=10.0_wvalues.txt")
#gs_saved = readdlm("../../plotting/data/N=10.0_gvalues.txt")
#gs = gs_saved / (10-1)

#read param values
#data_dir = "./phasediag_data_N10/"
data_dir = "../../data/phasediag_data_N10_zoomed/"
#mkpath(data_dir)
figs_dir = "../../paper_figs/"
#mkpath(figs_dir)
saved_files = readdir(data_dir)
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[1:end,1] #cut some off the end 
gs = params[1:end,2]

diffeq = (alg = Tsit5(),dt=1e-2,adaptive=false,maxiters=1e15)
diffeq = (alg = RK4(),dt=1e-2,adaptive=false,maxiters=1e15)

N = 10
#=
ws = range(0.5,0.85;length=50)
g_rel = 5
gs = range(1.2,3.0;length=50) .* (N-1)
=#
#=
#try original equations
u0 = generate_initial_condition(N)
rc=9.0;rp=9.0;K=1.0;c=2.4;α=0.005
p = [rc,rp,K,c,α]
ds = create_original_ds(u0,p)
=#

#critical params
wc1 = (N-1)/gs[1]
wc2 = (N-1)/N

#-------------------------create ds---------------------------------
#generate standard init cond
#u0 = rand(2N)
u0 = generate_initial_condition(N)
ds = create_ds(u0,[5,0.5,N];diffeq=diffeq) #no noise

σ = 1e-4
p_noisy = [5, 0.5, N, σ]
ds_noise = create_ds_noisy(u0,p_noisy;diffeq=diffeq) #noise 
#figure out how to calc lyap exp for noisy ds 
#d0_lower,d0_upper should be adjsted 



#---------------------time params----------------------
T_fp = 100
Ttr_fp = 100
Δt_fp = 0.1
fp_settings = [T_fp,Ttr_fp,Δt_fp]
Ttr = 2000
Ttrs = [Ttr]
T = 5000
lyap_vals = []

#calc lyapunov exponents
@time for (i,Ttr) in enumerate(Ttrs)
    @show Ttr
    @time Λ = lyapunov_diagram(ds_noise,gs,ws; T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
        d0=1e-3, ut=1e-1, lt=1e-4,save=false,realtime_plotting=false,parallelize=false)

    #@time Λ = lyapunov_diagram(ds,gs,ws; T=T, Ttr=Ttr, Δt=0.1, fp_settings=fp_settings, 
    #d0=1e-9,save=false,realtime_plotting=false,parallelize=false)
    push!(lyap_vals,Λ)
end


#=
#check lyap
ds = create_ds(u0, [18, 0.6, N]; diffeq=diffeq) #no noise
traj, = trajectory(ds, 0; Ttr=1000)
u0 = traj[end] #rand(2N) .* 0.1 
d0 = 1e-6
pds = ParallelDynamicalSystem(ds, [u0 , u0 .+ rand(2N) .* d0/sqrt(2N) ])
nsteps = 15000
d = zeros(nsteps)
for i in 1:nsteps
    step!(pds)
    states = current_states(pds)
    d[i] = norm(states[1] .- states[2])
end
pl = plot(d,lc=:gray10,yscale=:log10,
    label="no noise")

#check noisy lyap
#noise
ds_noise = create_ds_noisy(u0, [18, 0.6, N, 1e-4]; diffeq=diffeq) #noise 
u0 = traj[end] #rand(2N) .* 0.1
d0 = 1e-3
pds = ParallelDynamicalSystem(ds_noise, [u0, u0 .+ rand(2N) .* d0 / sqrt(2N)])
nsteps = 15000
d = zeros(nsteps)
for i in 1:nsteps
    step!(pds)
    states = current_states(pds)
    d[i] = norm(states[1] .- states[2])
end
plot!(pl, d, lc=:red, yscale=:log10, 
    label=L"σ=1e-4")
=#


savedir="../../data/phasediag_data_N10_noise_1e-4/"
mkpath(savedir)
writedlm(savedir*"phasediag_params.txt",hcat(ws,gs))
writedlm(savedir*"lyap_values.txt",lyap_vals)



#--------------read lyaps--------------------
data_dir = "../../data/phasediag_data_N10_no_noise/"
lyap_vals = reshape(readdlm(data_dir*"lyap_values.txt"),(50,50))
params = readdlm(data_dir*"phasediag_params.txt")
ws = params[:,1]
gs = params[:,2]
N=10


pl1 = plot_phase_diagram(lyap_vals[1],ws,gs,N;HEATMAP_DEFAULTS...,
    ylabel=L"\frac{1}{\epsilon(N-1)}",
    right_margin=6Plots.mm,
    left_margin=6Plots.mm,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_medium,
    xguidefontsize=guide_fontsize_large)
pl2 = plot_phase_diagram(lyap_vals[2],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
pl3 = plot_phase_diagram(lyap_vals[3],ws_saved[:],gs_N[:],N;HEATMAP_DEFAULTS...)
#savefig(pl1,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])_noise_1e-2"*".png")
savefig(pl1,figs_dir*"N10_phasediagram_Ttr_max_$(Ttrs[1])"*".png")

#savefig(pl2,"../../N5_phasediagram_Ttr_max_$(Ttrs[2])"*".png")
#savefig(pl3,"../../N5_phasediagram_Ttr_max_$(Ttrs[3])"*".png")


plot(ws,lyap_vals,ylims=(-0.2,Inf),labels = ["Ttrs=$(Ttrs[1])" "Ttrs=$(Ttrs[2])" "Ttrs=$(Ttrs[3])"])
