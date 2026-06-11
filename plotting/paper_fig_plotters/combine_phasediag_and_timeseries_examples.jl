using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

include("plot_phase_diag.jl")
include("plot_delay_timeseries_example.jl")

#save timeseries separately
pl_all_ts = plot(plots_ts...,layout=(1,4),size=(700,150), dpi = 300)
savefig(pl_all_ts,"../../paper_figs/fig1_N10_timeseries_examples.png")

#save timeseries separately
savefig(pl_phase_diag,"../../paper_figs/fig1_N10_phase_diag_Ttr1e4.png")


#combine phasediagram with time series
l = @layout [a{0.8h}; a{0.2w} a{0.3w} a{0.3w} a{0.2w} ]
pl_fig1 = plot(pl_phase_diag,plots_ts...,layout=l,size=(700,600))
savefig(pl_fig1,"../../paper_figs/fig1_phasediag_N10_timeseries_examples.png")