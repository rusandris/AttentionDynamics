using Revise
using AttentionDynamics
using OrdinaryDiffEq
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)


N = 300
K = 1.0
c = 2.4
r = 10.0
α0 = 0.005
#------------------------------solver settings----------------------------------
dtmax = 1e-2 #5*1e-2
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=dtmax, adaptive=true, maxiters=1e8)
#diffeq = (alg=Rodas5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e8)
#diffeq = (alg=Tsit5(), dt=1e-3, dtmax=5*1e-3, adaptive=true, maxiters=1e8)

T = 1000
Ttr = 2000
Δt = 0.1
nbins = 20

#create classic ds
x0 = rand(N) * 0.01
u0 = zeros(2N)
u0[1:2:end] .= x0

ds = create_original_ds(u0, (r, r, K, c, α0); diffeq=diffeq)

traj_classic, times_classic = trajectory(ds, T; Ttr=Ttr, Δt=Δt)

trajm_classic = Matrix(traj_classic)[:, 1:2:end]
pl_classic = plot(times_classic, trajm_classic, xlims=(Ttr + T - 50, Ttr + T), legend=false, title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")
plot(times_classic[1:500], trajm_classic[1:500, :], legend=false, yscale=:log10, ylimits=(1e-4, 1e1))

#integrate activity segment to segment
chart_time = 0.5
k = 10
int_act, time_segments = integrate_activity(trajm_classic, collect(times_classic), Int(chart_time / Δt))
lts_classic = get_lifetimes_integ(int_act, k)
pl_lt_classic = histogram(lts_classic,xlabel=L"\tau_{active} \cdot \Delta t", title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")

thresh = 0.01
lengths_act_classic = activity_periods(trajm_classic, thresh)
pl_lt_classic = histogram(lts_classic, bins=nbins, xlabel=L"\tau", title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")
pl_ap_classic = histogram(lengths_act_classic, bins=nbins, xlabel=L"\tau_{active} \cdot \Delta t", title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")


#------------------try varied delay system---------------------
#setup varied delays
αs = rand(N) .* 0.01 .+ 1e-3
ds_vd = create_original_varied_delay(u0, (r, r, K, c, αs); diffeq=diffeq)

#test trajectory

traj_vd, times_vd = trajectory(ds_vd, T; Ttr=Ttr, Δt=Δt)

trajm_vd = Matrix(traj_vd)[:,1:2:end]
pl_vd = plot(times_vd, trajm_vd, legend=false,
    title=L"r_p=%$(r),r_c=%$(r),\alpha \in (0.001,0.01)", xlims=(Ttr + T - 50, Ttr + T))
#pl_vd = plot(times_vd[1:100], trajm_vd[1:100,:], legend=false,yscale=:log10,ylimits=(1e-4,1e1))

#integrate activity segment to segment
chart_time = 0.5
k = 10
int_act, time_segments = integrate_activity(trajm_vd, collect(times_vd), Int(chart_time / Δt))
lts_vd = get_lifetimes_integ(int_act, k)
pl_lt_vd = histogram(lts_vd, bins=nbins, xlabel=L"\tau", title=L"r_p=%$(r),r_c=%$(r),\alpha \in (0.001,0.01)")

lengths_act_vd = activity_periods(trajm_vd, thresh)
pl_ap_vd = histogram(lengths_act_vd, bins=nbins, xlabel=L"\tau_{active} \cdot \Delta t", title=L"r_p=%$(r),r_c=%$(r),\alpha \in (0.001,0.01)")

#------------------try varied rate system---------------------
#rps = rand(6.0:1e-2:20.0, N)
rps = randn(N) .* 1 .+ 10
ds_vr = create_original_varied_rate(u0, (rps, r, K, c, α0); diffeq=diffeq)

traj_vr, times_vr = trajectory(ds_vr, T; Ttr=Ttr, Δt=Δt)

trajm_vr = Matrix(traj_vr)[:, 1:2:end]
pl_vr = plot(times_vr, trajm_vr, legend=false,
    title=L"r_p \in N(10,1) ,r_c=%$(r),\alpha = %$(α0)", xlims=(Ttr+T-50, Ttr+T))

chart_time = 0.5
k = 10
int_act, time_segments = integrate_activity(trajm_vr, collect(times_vr), Int(chart_time / Δt))
lts_vr = get_lifetimes_integ(int_act, k)

lengths_act_vr = activity_periods(trajm_vr .* Δt, thresh)
pl_ap_vr = histogram(lengths_act_vr, bins=nbins, xlabel=L"\tau_{active} \cdot \Delta t", title=L"r_p \in N(10,1) ,r_c=%$(r),\alpha = %$(α0)")
pl_lt_vr = histogram(lts_vr, bins=nbins, xlabel=L"\tau_{active}", title=L"r_p \in N(10,1) ,r_c=%$(r),\alpha = %$(α0)")

#------------------try asymmetric coupling---------------------

#setup interaction matrix
C = randn(N, N) .* 0.5 .+ 2.0
for i in 1:N
    C[i, i] = 0.0
    null_idxs = rand(1:N, rand(50:250))
    C[i, null_idxs] .= 0.0
    #anti_idxs = rand(1:N, rand(1:50))
    #C[i, anti_idxs] .= C[i, anti_idxs] .* (-1.0)
end

ds_as = create_original_asymmetric_interaction(u0, (r, r, K, C, α0); diffeq=diffeq)

traj_as, times_as = trajectory(ds_as, 50; Ttr=0, Δt=Δt)
trajm_as = Matrix(traj_as)[:, 1:2:end]

rs = [6.0,8.0,10.0]
p_sets_r = [[rs[i], rs[i]] for i in 1:length(rs)]
nr_ics=10
savedir_traj = "../../data/trajs_precomp_asymmetric/"
mkpath(savedir_traj)

#test precompute traj
precompute_trajectories(ds_as, nr_ics, p_sets_r, [1,2]; Δt=Δt_traj,
    T_traj=T, Ttr=Ttr, savedir=savedir_traj, p_names=["r"],
    ic_counter=ic_counter,
    filename="asymmetric_coupling_N$(N)_rp$(rs[2])_rp$(rs[2])_K$(K)_alpha$(α)_T$(T)_Ttr$(Ttr)_deltat$(Δt_traj).dat")


chart_time = 0.5
k = 10
int_act, time_segments = integrate_activity(trajm_as, collect(times_as), Int(chart_time / Δt))
lts_as = get_lifetimes_integ(int_act, k)

thresh = 0.01

lengths_act_as = activity_periods(trajm_as, thresh)
pl_lt_as = histogram(lts_as, bins=nbins, xlabel=L"\tau", title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")
pl_ap_as = histogram(lengths_act_as, bins=nbins, xlabel=L"\tau_{active} \cdot \Delta t", title=L"r_p=%$(r),r_c=%$(r),K=%$(K),c=%$(c),\alpha=%$(α0)")


#------------------------------------make plots-------------------------------------------

pl_lt = plot(pl_lt_classic,pl_lt_vd,pl_lt_vr,layout=(3,1),size=(400,800),dpi=300)
pl_ap = plot(pl_ap_classic, pl_ap_vd, pl_ap_vr, layout=(3, 1), size=(400, 800), dpi=300)
pl_ts = plot(pl_classic, pl_vd, pl_vr, layout=(3, 1), size=(400, 800), dpi=300)

savefig(pl_lt, "../../paper_figs/" * "test_lifetimes_histograms_varied_systems_alpha_var_rp_var_T$(T)_Ttr$(Ttr).png")
savefig(pl_ap, "../../paper_figs/" * "test_active_period_histograms_varied_systems_alpha_var_rp_var_T$(T)_Ttr$(Ttr).png")
savefig(pl_ts, "../../paper_figs/" * "test_lifetimes_varied_systems_alpha_var_rp_var_T$(T).png")