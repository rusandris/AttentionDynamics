using Plots
using LaTeXStrings
using DelimitedFiles
using DynamicalSystems
cd(@__DIR__)

data_dir = "../../data/data_escape_times/"
escape_time_files = readdir(data_dir)
plots_trs = []
ylims_trs = (0, 1800)


avg_escape_time_file = escape_time_files[findfirst(f -> occursin("escape_time", f), escape_time_files)]
avg_escape_time_data = readdlm(data_dir * avg_escape_time_file)
w_vals = readdlm(data_dir * "w_vals.txt")
Δws = avg_escape_time_data[:,1]
escape_times = avg_escape_time_data[:, 2]

#plot escape time scaling
pl_esc = plot(Δws, escape_times;
    xlabel=L"\Delta w", ylabel=L" \langle \tau_\mathrm{tr} \rangle", 
    xformatter=:auto,guidefontsize=20, widen=true,st=:scatter,
    ms=3, legend_columns=2,legend=false,yscale=:log10,xscale=:log10,
    mc=:gray20, markerstrokewidth=0.0, xlims=(0.17, 0.25), 
    xticks=[10^(-0.8), 10^(-0.7), 10^(-0.6)],
    tickfontsize=10)

#plot transient decay examples 
#from individual parameters 
ensemble_decay1 = readdlm(data_dir * "lt_ts_0.5")
ensemble_decay2 = readdlm(data_dir * "lt_ts_0.49")
ensemble_decay3 = readdlm(data_dir * "lt_ts_0.48")

plot!(inset=(1, bbox(0.55, 0.07, 0.4, 0.4)), subplot=2)
plot!(pl_esc[2], ensemble_decay1[1:end-1, 1], ensemble_decay1[1:end-1, 2],lc=:blue,yscale=:log10,legend=false,
    ylabel=L"N_e(t)", xlabel=L"t", guidefontsize=15,yticks=[1e1,1e3,1e5], xticks=[0, 5000])
plot!(pl_esc[2], ensemble_decay2[1:end-1, 1], ensemble_decay2[1:end-1, 2], lc=:orange, yscale=:log10)
plot!(pl_esc[2], ensemble_decay3[1:end-1, 1], ensemble_decay3[1:end-1, 2], lc=:green, yscale=:log10)

#replot example points with colors 
plot!(pl_esc[1], [w_vals[41] - 0.3], [escape_times[41]],
    st=:scatter,mc=:blue,ms=8,markershape=:diamond,
    yscale=:log10, xscale=:log10)

plot!(pl_esc[1], [w_vals[21] - 0.3], [escape_times[21]],
    st=:scatter, mc=:orange, ms=6, markershape=:rect,
    yscale=:log10, xscale=:log10)

plot!(pl_esc[1], [w_vals[1] - 0.3], [escape_times[1]],
    st=:scatter, mc=:green, ms=6, markershape=:utriangle,
    yscale=:log10, xscale=:log10)

#add fit line with offset
b, a = linreg(log10.(Δws), log10.(escape_times))
plot!(pl_esc[1], Δws, 10 .^ (b-0.1) .* Δws .^a,ls=:dash,lc=:red)