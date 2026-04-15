using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("defaults.jl")
logocolors = Colors.JULIA_LOGO_COLORS
figs_dir = "./poster_figs/"

#read attractor finding results
data_dir = "../../data/attractor_finding/nr_attractors_cumulative_avg/"
nr_atrs_avgs = readdlm(data_dir * "found_attractors_nr_cumulative_avg_N3N4N5.txt")
#method_params = readdlm("../../data/attractor_finding/attractor_finding_method_params.txt")
#Ns = Int.(method_params[1:3])
Ns = [3,4,5]

#--------------------------------plot attractor finding----------------------------
lcolors = [logocolors.red,logocolors.green,logocolors.purple]
markershapes = [:circle,:dtriangle,:utriangle]
pl_atr = plot(;ylabel="# attractors",xlabel="# initial conditions",
    legend=:bottomright,legend_columns=-1,DEFAULTS...,legendfontsize=legend_fontsize+5)

for (i,N) in enumerate(Ns)

    nr_atrs_avg = nr_atrs_avgs[:,i]
    #------------------plot----------------------------
    nr_atr_max = factorial(Ns[i]-1)
    annotate!(35, nr_atr_max*1.2, text(L"(N-1)! = %$nr_atr_max", :gray10, :center, 28))
    #plot!(nr_atrs_avg ./ length(resu   lts),lw=4,lc=lcolors[i],label=L"N=%$(Ns[i])")
    plot!(nr_atrs_avg[1:3:end], ms=8, st=:scatter, lc=:gray10,mc=:gray10,  
        markershape=markershapes[i], label=L"N=%$(Ns[i])")
    hline!([nr_atr_max],lw=3,ls=:dash,lc=:gray10,label="")

end
#yticks = logrange(1,40;length=4)
plot!(pl_atr,ylims=(1,40),xlims=(1,50),yticks=[1,10],yaxis=:log,tickfontsize=tick_fontsize_large+5,guidefontsize=guide_fontsize_large+5)
plot!(size=(800,700))


savefig(pl_atr,figs_dir*"system_multistability_nr_attractors.png")
#=
data_dir = "./attractor_finding/" 
mkpath(data_dir)
writedlm(data_dir*"attractor_finding_Ns345.txt",stack(nr_atrs_avgs))
writedlm(data_dir*"attractor_finding_method_params.txt",[Ns...,g_rel,w,cons_recurrences,locate_steps,max_ensemble,nr_runs])
=#