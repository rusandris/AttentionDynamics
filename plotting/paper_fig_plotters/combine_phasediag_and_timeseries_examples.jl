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
l = @layout [a{0.8w} grid(4,1)]
fig1 = plot(pl_phase_diag,plots_ts...,layout=l,size=(1000,500))