using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

include("plot_phase_diag_N10.jl")
include("plot_delay_timeseries_example.jl")

#save timeseries separately
l_ts = @layout [a{0.2w} a{0.3w} a{0.3w} a{0.2w} ]
pl_all_ts = plot(plots_ts...,layout = l_ts,size=(700,150), dpi = 300,top_margin=3Plots.mm)

#add fig level annotation
#annotate!(pl_all_ts[1], 0.0, 5, text(L"(a)", :gray10, :center, annot_size + 5))

savefig(pl_all_ts,"../../paper_figs/fig1_N10_timeseries_examples.png")

#save timeseries separately
savefig(pl_phase_diag,"../../paper_figs/fig1_N10_phase_diag_Ttr1e4.png")

#add fig level annotation 
annotate!(pl_phase_diag[end-1],-0.17, 10.9, text(L"(b)", :gray10, :center, annot_size))
annotate!(pl_all_ts[1], 0.0, 7.5, text(L"(a)", :gray10, :center, annot_size))

#combine phasediagram with time series

l = @layout [a{0.2h};b{0.8h}]
pl_fig1 = plot(pl_all_ts,pl_phase_diag,layout=l,size=(700,600))


savefig(pl_fig1,"../../paper_figs/fig1_phasediag_N10_timeseries_examples.png")