using Pkg
#cd("Documents/stn_research/STNResearch/fibrillation/")
cd(@__DIR__)
Pkg.activate(".")
using DelimitedFiles
using CSV,DataFrames
using Plots,LaTeXStrings
#plotlyjs()
gr()


#-------------------------------LTAF Database------------------------------
data = CSV.read("Long_Term_AF_Database/preAF_fib/sig_06.txt",DataFrame;types=[Float64,Float64,Float64],header=1,skipto=3)
rrs = readdlm("Long_Term_AF_Database/preAF_fib/rr_06.txt")

#-------------------------------NSR Database-------------------------------
#=
data = CSV.read("NSR_Database/NSR/sig_16265.txt",DataFrame;types=[Float64,Float64,Float64],header=1,skipto=3)
rrs = readdlm("NSR_Database/NSR/rr_16265.txt")
=#

times_elapsed = data[:,1] .- data[1,1]
rrs_start_times = Float64.(rrs[:,1])
rrs_start_times_elapsed = rrs_start_times .- rrs_start_times[1]

pl_sig = plot(times_elapsed[1:1000,1],data[1:1000,2],xlims=[1,5],
    xlabel=L"t [s]",ylabel=L"ECG",lc=:gray10,dpi=300,
    label="ECG signal 1",legend=false,guidefontsize=18, legendfontsize=12)

vline!(pl_sig,rrs_start_times_elapsed,lc=:red,ls=:dot,lw=1,
    ylims=[-1,4.5],label="RR annot")

#savefig(pl_sig,"sig_sample16_preAF_example.png")
savefig(pl_sig,"sig_sample06_LTAF_example.png")

rr_intervals = rrs[2:end,2]
pl_RR = plot(rr_intervals[:],
    markershape=:circle, ylims=(0.3,1.0),
    lc=:gray10, mc=:gray10, 
    yguidefontsize=18, legendfontsize=12, xguidefontsize=12,
    xlabel="heartbeat #",ylabel=L"\Delta t_{RR}[s]",legend=false,
    xlims=[2500,2700],dpi=300)#markershape=:circle,ms=2)

pl_RR_resized = plot(pl_RR,size=(800,300),left_margin=5Plots.mm,
    bottom_margin=5Plots.mm, dpi=300)


#savefig(pl_RR,"rr_sample16_preAF_example.png")
savefig(pl_RR_resized,"rr_sample06_LTAF_example.png")
