using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

include("defaults.jl")
include("../plot_functions/linreg.jl")

datadir = "../../data/scaling_data_N2/"
scaling_data = readdlm(datadir*"scalingdata_N2_g5.0.txt")

N = 2
g = 5.0

wc1 = (N-1)/g
wc2 = (N-1)/N

hc_params = scaling_data[:,1]
periods = scaling_data[:,2]
hopf_params = scaling_data[:,3]
amplitudes = scaling_data[:,4]

Δw_hopf = hopf_params .- wc1
Δw_clinic = wc2 .- hc_params

#-----------------------------hopf---------------------------------
ms = 2
pl1 = plot(Δw_hopf, amplitudes; DEFAULTS..., st=:scatter, ms=ms, 
    xlabel=L"\Delta w", ylabel=L"A(x_1)", legend=true, label="",
    mc=:gray10, tickfontsize=14, xticks = [0.0,0.05,0.1],
    yticks=[0.0, 0.5, 1.0,1.5],ylims=(0.0,1.5))
    #yscale=:log10,xscale=:log10)

a1,b1 = my_linreg(log10.(Δw_hopf),log10.(amplitudes))
@show a1,b1

#plot!(Δw1,(x -> 10^b1*x^a1).(Δw1);ls=:dash,lc=:red)
#plot!(yaxis=:log,xaxis=:log)

pl_hopf = plot!(Δw_hopf,(x -> 10^b1*sqrt(x)).(Δw_hopf);
    ls=:dash,lc=:red,label="",
    legendfontsize=DEFAULTS.xguidefontsize)
    #label=L"0.6 \cdot \Delta w^{1/2}")

#----------------------------homoclinic----------------------------------

a2,b2 = my_linreg(log10.(reverse(Δw_clinic)),log10.(reverse(periods)))
@show a2,b2
a2r = round(a2;digits=2)

pl2_lin = plot(reverse(Δw_clinic),reverse(periods);DEFAULTS...,st=:scatter,
    ms=ms, xlabel=L"\Delta w", ylabel=L"T", legend=true, 
    label="", mc=:gray10,ylims=(0,Inf))
plot!(Δw_clinic,(x -> 10^b2*x^a2).(Δw_clinic);ls=:dash,lc=:red,
    legendfontsize=DEFAULTS.xguidefontsize,mc=:gray10,label="")
    #label=L"0.035 \cdot \Delta w^{-0.84}")
#pl2 = plot!(yaxis=:log,xaxis=:log,xticks=[10^(-2.5),10^(-2),10^(-1.5)])

pl_clinic = plot(log10.(reverse(Δw_clinic)), log10.(reverse(periods)); DEFAULTS..., st=:scatter,
    ms=ms, xlabel=L"\log(\Delta w)", ylabel=L"\log(T)",
    legend=true, label="", mc=:gray10,ylims=(0,Inf))

plot!(log10.(Δw_clinic), log10.((x -> 10^b2 * x^a2).(Δw_clinic)); 
    ls=:dash, lc=:red,legendfontsize=DEFAULTS.xguidefontsize, 
    mc=:gray10, label="",ticks=(-3:3,[L"10^{%$(i)}" for i in -3:3]),
    tickfontsize=18)

savefig(pl_hopf, "../../paper_figs/hopf_amplitude_scaling_N2_g5.png")
savefig(pl_clinic, "../../paper_figs/homoclinic_period_scaling_N2_g5.png")


