using Revise
using ComplexityMeasures
using BifurcationKit
using AttentionDynamics
using OrdinaryDiffEq
using StateTransitionNetworks
using DynamicalSystems
using StatsBase
using Plots, LaTeXStrings
using LinearAlgebra
using DelimitedFiles
cd(@__DIR__)

#------------------IDEA---------------------
#test multidim OP after critical point 
#what happens to the symbolic dynamics of w is changed?
#PE,WPE
#distribution / average lifetime of in the top k list 
#changing a clearly interpretable parameter: α,c 
#corresponds to moving along a curve in g-w space
#which parameter change leads to driving the system more towards criticality?
#how does properties of the system change when it gets near criticality?


#TODO:
#try exponential/statistically significant binning 
#change charting time (3plot) + write time sampling on axis
#charting time 1/10th of average peak width
#add fits
#show inequality on fig
#separate lifetimes from plotting and save 

include("chart_functions.jl")

gauss(x,μ,σ) = 1/(σ*sqrt(2π)) * exp(-(x - μ)^2/(2σ^2))

nbins = 50
figsize = (600,800)
legendfontsize=10
savedir = "../../paper_figs/lifetime_distributions/test6/"
mkpath(savedir)
data_dir_traj_root = "../../data/trajs_precomp_test6/"
data_dir_lt_root = "../../data/lifetimes_precomp_test6/"


#----------------------------change r parameter-------------------------
# try rp = rc scenario
# compare with rp \neq rc (mismatch between production and consumption rate)
# can the model produce different distributions?
# try different k 
# try different Δt

Δt = 1e-1
K = 1.0;
c = 2.4;
α = 0.005;

#READ trajectories
data_dir = data_dir_traj_root*"rp=rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)[2:end]
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end
#read params 
rs = vec(readdlm(data_dir*files_all[1]))

#get lifetimes with different k, charting times 
ks = [10,50,100]
Δts_chart = [1,10,30]

#precompute lifetimes 
data_dir = data_dir_lt_root * "rp=rc/"
mkpath(data_dir)
for (k_idx, k) in enumerate(ks)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, k, Δts_chart[1])

    #save lifetimes 
    for (i, lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rp=rc=$(rs[i])_k_$(k)_Δt$(Δts_chart[1]).txt", lts)
    end

end

data_dir = data_dir_lt_root * "rp=rc/"
mkpath(data_dir)
for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)

    #save lifetimes 
    for (i,lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rp=rc=$(rs[i])_k_$(ks[1])_Δt$(Δt).txt", lts)
    end
    
end


# compare with rp \neq rc (mismatch between production and consumption rate)

#----------------------------------rp----------------------------------------

#READ trajectories
data_dir = data_dir_traj_root * "/rp/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt",f),files_all)[2:end]
trajs = Matrix{Float64}[]
for i in 1:length(traj_files) 
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir*traj_files[i],Float64)
    push!(trajs,traj)
end
#read rs
rs = vec(readdlm(data_dir * files_all[1]))
rc = parse(Float64,files_all[1][end-5:end-4])


data_dir = data_dir_lt_root * "rp/"
mkpath(data_dir)
for (k_idx, k) in enumerate(ks)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, k, Δts_chart[1])

    #save lifetimes 
    for (i, lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rp$(rs[i])_rc$(rc)_k_$(k)_Δt$(Δts_chart[1]).txt", lts)
    end

end

data_dir = data_dir_lt_root * "rp/"
mkpath(data_dir)
for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)

    #save lifetimes 
    for (i, lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rp$(rs[i])_rc$(rc)_k_$(ks[1])_Δt$(Δt).txt", lts)
    end

end


#WATCH OUT FOR NUMERICAL METHOD ERRORS
# compare with rp \neq rc (mismatch between production and consumption rate)
#ONLY CHANGE rc
#precompute trajectories 


#----------------------------rc----------------------------------

#READ trajectories
data_dir = data_dir_traj_root * "/rc/"
files_all = readdir(data_dir)
traj_files = filter(f -> occursin(".txt", f), files_all)[2:end]
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir * traj_files[i], Float64)
    push!(trajs, traj)
end

#read rs
rs = vec(readdlm(data_dir * files_all[1]))
rp = parse(Float64, files_all[1][end-4:end-3])

data_dir = data_dir_lt_root * "rc/"
mkpath(data_dir)
for (k_idx, k) in enumerate(ks)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, k, Δts_chart[1])

    #save lifetimes 
    for (i, lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rc$(rs[i])_rp$(rp)_k_$(k)_Δt$(Δts_chart[1]).txt", lts)
    end

end

data_dir = data_dir_lt_root * "rc/"
mkpath(data_dir)
for (Δt_idx, Δt) in enumerate(Δts_chart)
    lifetimes_x, lifetimes_y = get_lifetimes(trajs, ks[1], Δt)

    #save lifetimes 
    for (i, lts) in enumerate(lifetimes_x)
        writedlm(data_dir * "lifetimes_x_rc$(rs[i])_rp$(rc)_k_$(ks[1])_Δt$(Δt).txt", lts)
    end

end
