using Plots,LaTeXStrings
using DelimitedFiles

include("defaults.jl")
include("plot_functions/linreg.jl")

scaling_data = readdlm("data/scalingdata_N2_g5.0.txt")

N = 2
g = 5.0

wc1 = (N-1)/g
wc2 = (N-1)/N

hc_params = scaling_data[:,1]
periods = scaling_data[:,2]
hopf_params = scaling_data[:,3]
amplitudes = scaling_data[:,4]

Δw1 = hopf_params[2:end] .- wc1
Δw2 = wc2 .- hc_params

pl1 = plot(Δw1,amplitudes;DEFAULTS...,st=:scatter,ms=3,xlabel=L"\Delta w",ylabel=L"A(x_1)")
a1,b1 = my_linreg(log10.(Δw1),log10.(amplitudes[2:end]))
@show a1,b1

#plot!(Δw1,(x -> 10^b1*x^a1).(Δw1);ls=:dash,lc=:red)
#plot!(yaxis=:log,xaxis=:log)

pl1 = plot!(Δw1,(x -> 10^b1*sqrt(x)).(Δw1);ls=:dash,lc=:red)
savefig(pl1,"hopf_amplitude_scaling_N2_g5.png")


a2,b2 = my_linreg(log10.(reverse(Δw2)),log10.(reverse(periods)))
@show a2,b2

pl2 = plot(reverse(Δw2),reverse(periods);DEFAULTS...,st=:scatter,ms=3,xlabel=L"\Delta w",ylabel="T")
pl2 = plot!(Δw2,(x -> 10^b2*x^a2).(Δw2);ls=:dash,lc=:red)
pl2 = plot!(yaxis=:log,xaxis=:log)

savefig(pl2,"homoclinic_period_scaling_N2_g5.png")


