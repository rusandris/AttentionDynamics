using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
using DataFrames, CSV
using StatsBase
using CurveFit  
cd(@__DIR__)

include("plot_phase_diag_N.jl")
include("plot_lyap_curves_largesys.jl")

pl_fig1 = plot(pl_phasediag_N,pl_lyaps,layout=(1,2),size=(1800,600),
    guidefontsize=25,tickfontsize=15,left_margin=18Plots.mm,bottom_margin=10Plots.mm,
    legendfontsize=18,titlefontsize=22)