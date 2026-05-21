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

#combine phasediagram with time series
l = @layout [grid(1,4); a{0.6h}]
fig1 = plot(plots_ts...,pl_phase_diag,layout=l,size=(1000,800))