
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

savedir = "../../paper_figs/lifetime_distributions/test10/"
mkpath(savedir)
data_dir_traj_root = "../../data/trajs_precomp_test_c/"
data_dir_lt_root = "../../data/lifetimes_precomp_test_c/"


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


#----------------------------c----------------------------------
Δt = 5*1e-2
#READ trajectories
data_dir_traj_root = "../../data/trajs_precomp_test_c/"
data_dir_traj = data_dir_traj_root * "c_long/"
files_all = readdir(data_dir_traj)
traj_files = filter(f -> occursin(".txt", f), files_all)[2:end]

#=
#test one-by-one
traj = readdlm(data_dir_traj * traj_files[1], Float64)
#instant
lifetimes_c, _ = get_lifetimes(traj, 20, 5)

pl_distrib_c = plot_distrib(lifetimes_c, 5;
    color=:gray10, nbins=40, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")


#integrated
ts = [0.0:Δt:10000;]
L_ints, time_segments = integrate_activity(traj, ts, 5)
lts_integ = get_lifetimes_integ(L_ints, 20)

pl_distrib_c_integ = plot_distrib(lts_integ, 5;
    color=:gray10, nbins=50, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")

=#


#=
#read trajs together
trajs = Matrix{Float64}[]
for i in 1:length(traj_files)
    println("Reading file " * traj_files[i] * "...")
    traj = readdlm(data_dir_traj * traj_files[i], Float64)
    push!(trajs, traj)
end
=#

#read cs
cs = vec(readdlm(data_dir_traj * files_all[1]))
r = parse(Float64, files_all[2][10:12])
Δt = parse(Float64, files_all[2][end-6:end-4])


data_dir = data_dir_lt_root * "c_long/"
mkpath(data_dir)
for (k_idx, k) in enumerate(ks)
    @show k
    #read each traj at a time 
    for i in 1:length(traj_files)
        @show i
        traj = readdlm(data_dir_traj * traj_files[i], Float64)
        lifetimes_x, _ = get_lifetimes(traj, k, Δts_chart[1])
        writedlm(data_dir * "lifetimes_x_rc$(r)_rp$(r)_k_$(k)_Δt_$(Δt)_chart_time_$(Δts_chart[1])_c$(cs[i]).txt", lifetimes_x)

    end
    
end

ks = [10,50,100]
Δts_chart = [1,5,10]
ts = [0.0:Δt:10000;]

k=20
data_dir = data_dir_lt_root * "c_long/"
mkpath(data_dir*"/instant/")
mkpath(data_dir*"/integrated/")
for (Δt_idx, time_chart) in enumerate(Δts_chart)
    @show time_chart
    #read each traj at a time 
    for i in 1:length(traj_files)
        @show i
        traj = readdlm(data_dir_traj * traj_files[i], Float64)

        #instant
        lifetimes_x, _ = get_lifetimes(traj, k, time_chart)
        writedlm(data_dir * "/instant/lifetimes_x_rc$(r)_rp$(r)_k_$(k)_Δt_$(Δt)_chart_time_$(time_chart)_c$(cs[i]).txt", lifetimes_x)

        L_ints, time_segments = integrate_activity(traj, ts, time_chart)
        lts_integ = get_lifetimes_integ(L_ints, k)
        writedlm(data_dir * "/integrated/lifetimes_x_rc$(r)_rp$(r)_k_$(k)_Δt_$(Δt)_chart_time_$(time_chart)_c$(cs[i]).txt", lts_integ)

    end

end

#------------------try varied delay system---------------------
data_dir_traj_root = "../../data/trajs_precomp_test_varied_delay/"
data_dir_traj = data_dir_traj_root * "1e-6_1e-2/"
files_all = readdir(data_dir_traj)
traj_files = filter(f -> occursin(".txt", f), files_all)
r = parse(Float64, files_all[1][10:12])


traj = readdlm(data_dir_traj * traj_files[1], Float64)
lifetimes_x1, _ = get_lifetimes(traj, 10, 1)


pl_distrib_x1 = plot_distrib(lifetimes_x1, 1;
    color=:gray10, nbins=50, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")

#try integrated variables
ts = [0.0:Δt:10000;]
L_ints, time_segments = integrate_activity(traj, ts, 10)
lts_integ_x1 = get_lifetimes_integ(L_ints, 20)

pl_distrib_x1_integ = plot_distrib(lts_integ_x1, 5;
    color=:gray10, nbins=40, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")


#------------------try varied rate system---------------------
data_dir_traj_root = "../../data/trajs_precomp_test_varied_rate/"
data_dir_traj = data_dir_traj_root * "6_20/"
files_all = readdir(data_dir_traj)
traj_files = filter(f -> occursin(".txt", f), files_all)
r = parse(Float64, files_all[1][10:12])


traj = readdlm(data_dir_traj * traj_files[2], Float64)
lifetimes_x2, _ = get_lifetimes(traj, 10, 5)


pl_distrib_x2 = plot_distrib(lifetimes_x2, 1;
    color=:gray10, nbins=30, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")

#------------------try varied rate, delay system---------------------
data_dir_traj_root = "../../data/trajs_precomp_test_varied_rate_delay/"
data_dir_traj = data_dir_traj_root
files_all = readdir(data_dir_traj)
traj_files = filter(f -> occursin(".txt", f), files_all)
r = parse(Float64, files_all[1][9:11])
Δt = parse(Float64, files_all[1][end-6:end-4])

traj = readdlm(data_dir_traj * traj_files[1], Float64)
lifetimes_x3, _ = get_lifetimes(traj, 20, 1)

pl_distrib_x3 = plot_distrib(lifetimes_x3, 1;
    color=:gray10, nbins=30, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")

#try integrated variables
ts = [0.0:Δt:10000;]
L_ints, time_segments = integrate_activity(traj, ts, 5)
lts_integ = get_lifetimes_integ(L_ints, 20)

pl_distrib_x3_integ = plot_distrib(lts_integ, 5;
    color=:gray10, nbins=50, mode=:none,
    scale=:log,
    title=L"k=10" * L",Δt=%$(1*Δt)",
    label_par=L"r_c = $(r), r_p = $(r)")
