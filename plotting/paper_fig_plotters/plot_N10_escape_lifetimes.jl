using Plots
using LaTeXStrings
using AttentionDynamics
using DelimitedFiles
using DynamicalSystems
using OrdinaryDiffEq
cd(@__DIR__)

#-------------------------read transient time files-------------------------------

data_dir_root = "../../data/transient_lifetimes_from_ensembles_N10_g45_nens1000_Tmax1e5_merged_recalc_end_Tmax1e6/"
data_dir = data_dir_root * "transient_lifetimes_from_ensembles_Tmax1e5_merged_recalc_end_Tmax1e6/"
escape_time_files = readdir(data_dir)

#--------------------------make escape number curves----------------------------
Δt_bin = 1.0
N0 = 1e3
@info "Making transient escape number curves..."
#loop over transient (escape) time files 
transient_curves = []
curve_slopes_interc = []
nr_infs = [] 
plots_slopes = []
for i in 1:length(escape_time_files)

    #readfile
    file = escape_time_files[i]
    transient_times = readdlm(data_dir*file)

    #count infs for validation
    cnt_infs = count(t -> isinf(t), transient_times)
    push!(nr_infs, cnt_infs)

    #make curve by binning
    Ne, ts = transient_decay_curve(transient_times; Δt=Δt_bin)
    push!(transient_curves, (Ne, ts))

    #adaptive line fitting interval 
    slope_measuring_interval = round(Int, 0.05 * length(ts)):round(Int, 0.5 * length(ts))
    @show slope_measuring_interval
    #linear regression to decay curve
    b, a = linreg(ts[slope_measuring_interval], log10.(Ne[slope_measuring_interval]))
    push!(curve_slopes_interc,(a,b))

    #check slopes by plotting
    pl = plot(ts[1:end-1], Ne[1:end-1], lc=:gray10, yscale=:log10)
    plot!(pl, ts, 10^(b) .* 10 .^ (ts .* a), ls=:dash, lc=:red, yscale=:log10)
    push!(plots_slopes,pl)
end
valid_runs = nr_infs .< 100
transient_curves_valid = transient_curves[valid_runs]
curve_slopes_interc_valid = curve_slopes_interc[valid_runs]
slopes = [curve_slopes_interc[i][1] for i in 1:length(transient_curves_valid)]
escape_times = [-1/s for s in slopes] #from fitting
 
#-----------------------------------plot------------------------------------


w_vals = readdlm(data_dir_root * "w_vals_transients_N10_g45_nens1000_Tmax1e5.txt")
w_vals_valid = w_vals[valid_runs]
w_NS = 0.329206866618409
Δws = w_vals_valid .- w_NS
Δlog_ws = log10.(w_vals_valid) .- log10(w_NS)


#plot escape time scaling
pl_esc = plot(Δws, escape_times;
    xlabel=L"\Delta w", ylabel=L" \langle \tau_\mathrm{tr} \rangle", 
    xformatter=:auto,guidefontsize=20, widen=true,st=:scatter,
    ms=3, legend_columns=2,legend=false,yscale=:log10,xscale=:log10,
    mc=:gray20, markerstrokewidth=0.0, xlims=(0.0005, 0.15),
    ylims=(1e2, 5*1e5),yticks = [1e2,1e3,1e4,1e5], 
    tickfontsize=10, xticks=[1e-4, 1e-3, 1e-2, 1e-1])


#=
#plot escape time scaling
pl_esc = plot(Δlog_ws, log10.(escape_times);
    xlabel=L"\log w - \log w_{NS}", ylabel=L" \log(\langle \tau_\mathrm{tr} \rangle)",
    xformatter=:auto, guidefontsize=20, widen=true, st=:scatter,
    ms=2, legend_columns=2, legend=false, ylims=(2,5),xlims=(-0.07,0.16),
    mc=:gray20, markerstrokewidth=0.0,tickfontsize=10)
=#

#-----------------------plot transient decay examples---------------------------
#from individual parameters 
Ne1,ts_e1 = transient_curves[25]
Ne2, ts_e2 = transient_curves[30]
Ne3, ts_e3 = transient_curves[35]

plot!(inset=(1, bbox(0.17, 0.45, 0.35, 0.35)), subplot=2)
plot!(pl_esc[2], ts_e1[1:end-1], Ne1[1:end-1], lc=:blue, yscale=:log10, legend=false,
    ylabel=L"N_s(t)", xlabel=L"t", guidefontsize=15,yticks=[1e1,1e3,1e5], xticks=[0, 8000])
plot!(pl_esc[2], ts_e2[1:end-1], Ne2[1:end-1], lc=:orange, yscale=:log10)
plot!(pl_esc[2], ts_e3[1:end-1], Ne3[1:end-1], lc=:green, yscale=:log10)

#replot example points with colors 
plot!(pl_esc[1], [Δws[25]], [escape_times[25]],
    st=:scatter,mc=:blue,ms=8,markershape=:diamond)

plot!(pl_esc[1], [Δws[30]], [escape_times[30]],
    st=:scatter, mc=:orange, ms=6, markershape=:rect)

plot!(pl_esc[1], [Δws[35]], [escape_times[35]],
    st=:scatter, mc=:green, ms=6, markershape=:utriangle)


#add fit line with offset
b, a = linreg(log10.(Δws[1:end-20]), log10.(escape_times[1:end-20]))
plot!(pl_esc[1], Δws[1:end-20], 10 .^ (b + 0.3) .* Δws[1:end-20] .^ a, ls=:dash, lc=:red)


#--------------------------plot timeseries transition point--------------------------
#inset 2
N=10
w_example_ts = w_vals[25] 

p = [45, w_example_ts, N]
u0 = rand(2N)

diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e8)
#u0 = rand(2N)
u0 = zeros(2N)
u0[1:2:end] = [0.1:0.1:1.0;] 
ds = create_ds(u0, p; diffeq=diffeq)

τ_tr = simple_transient_time(ds, u0; Tmax=4000, threshold=1.0,period_tol=1e-3,max_tol=1e-3)
traj, t = trajectory(ds, τ_tr+100, u0; Ttr=0, Δt=1e-2)

#plot time series 
trajm = Matrix(traj)[:,1:2:end]
plot!(pl_esc, inset=(1, bbox(0.56, 0.07, 0.45, 0.2)), subplot=3)
plot!(pl_esc[3], t, trajm, legend=false, lw=0.5, xlims=(τ_tr-120, τ_tr+40), ylims=(0, 40),
    ylabel="", xlabel=L"t", yticks=[0,20,40], guidefontsize=15)

#adjust limits
#plot!(pl_esc[3],xlims=(100,300))

#add annotation
annotate!(pl_esc[3], τ_tr - 155, 20, text(L"x_i", :gray10, :center, 18))
