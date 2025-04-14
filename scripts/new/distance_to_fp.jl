using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)

function dist_to_fixedpoint(u,fp)
    return norm(u - fp)
end


#set parameters
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

N = 6
#N = 3
#N = 5 
g = g_N * (N-1) 
wc1 = (N-1)/g
wc2 = (N-1)/N
w_test = 0.3 #wc1+ (wc2-wc1)/2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

ds = create_ds(rand(2N),[g,w_test,N];diffeq=diffeq)
traj,t = trajectory(ds,300;Ttr=1000)



function fixed_point(N,M,p) #use [x1,y2,x2,y2...] ordering
    fp = zeros(2N) #init fixed point state
    w = p[2]
    w_M = (w*M)/(N-1)
    w_N = (w*N)/(N-1)
    s = 1/(1 + w_M - w_N) #fixed point value (nonzero)
    fp[1:2:2M] .= s #set x vars to s
    fp[2:2:2M] .= s #set y vars to s
    return fp
end

#fp_coex = fixed_point(N,N,p)
#fp_dom = fixed_point(N,1,p)
#fp_M2 = fixed_point(N,2,p)

ws = [0.3:0.0025:0.9;]

λ_specs_coex = []
λ_specs_dom = []
λ_specs_M2 = []
λ_specs_M3 = []
λ_specs_M4 = []

dist_extremas_dom = []
dist_extremas_coex = []

for (i,w) in enumerate(ws)
    @show w
    #set parameters
    set_parameter!(ds,2,w)
    p = current_parameters(ds)
    @show p

    #recalc fixed points with new parameters
    fp_coex = fixed_point(N,N,p)
    fp_dom = fixed_point(N,1,p)
    fp_dom2 = reverse(fp_dom) #generate another valid dominance fp easily 

    traj,ts = trajectory(ds,200,initial_state(ds);Ttr=1000,Δt=1e-2)

    #distance between two dominance fixed points 
    d_dom_fp = norm(fp_dom-fp_dom2)

    @show traj[end]
    #calculate distances to fixed points along the trajectory
    dists_coex = [dist_to_fixedpoint(traj[i],fp_coex) for i in 1:length(traj)]

    dists_dom = [dist_to_fixedpoint(traj[i],fp_dom) for i in 1:length(traj)]    

    if w >= wc2 && !(isapprox(norm(dists_dom[end]),0,atol=1e-2))
        dists_dom = dists_dom .- d_dom_fp 
    end

    #calculate various fp spectra
    spec_coex = fixed_point_eigen_spectrum(N,N,p) #coexistence spectrum
    spec_dom = fixed_point_eigen_spectrum(N,1,p) #dominance spectrum
    spec_M2 = fixed_point_eigen_spectrum(N,2,p) #M=2 spectrum
    spec_M3 = fixed_point_eigen_spectrum(N,3,p) #M=3 spectrum
    spec_M4 = fixed_point_eigen_spectrum(N,4,p) #M=4 spectrum

    push!(λ_specs_coex,spec_coex)
    push!(λ_specs_dom,spec_dom)
    push!(λ_specs_M2,spec_M2)
    push!(λ_specs_M3,spec_M3)
    push!(λ_specs_M4,spec_M4)

    push!(dist_extremas_dom,extrema(dists_dom))
    push!(dist_extremas_coex,extrema(dists_coex))
end

#dom

dist_extremas_dom_mins = [dist_extremas_dom[i][1] for i in 1:length(dist_extremas_dom)] 
dist_extremas_dom_maxs = [dist_extremas_dom[i][2] for i in 1:length(dist_extremas_dom)]

#coex
dist_extremas_coex_mins = [dist_extremas_coex[i][1] for i in 1:length(dist_extremas_coex)]
dist_extremas_coex_maxs = [dist_extremas_coex[i][2] for i in 1:length(dist_extremas_coex)]

#extract real and imaginary parts 
realparts_coex = stack(real.(λ_specs_coex))'
imagparts_coex = stack(imag.(λ_specs_coex))'
realparts_dom = stack(real.(λ_specs_dom))'
imagparts_dom = stack(imag.(λ_specs_dom))'
realparts_M2 = stack(real.(λ_specs_M2))'
imagparts_M2 = stack(imag.(λ_specs_M2))'
realparts_M3 = stack(real.(λ_specs_M3))'
imagparts_M3 = stack(imag.(λ_specs_M3))'
realparts_M4 = stack(real.(λ_specs_M4))'
imagparts_M4 = stack(imag.(λ_specs_M4))'


ylims_spec = [-6.0,3.0]

#plot coexistence fixed point spectrum
pl = plot(ws,realparts_coex[:,[1,2]],ylabel="eigen spectrum",xlabel=L"w",legend=:outertopright,
    title="N=$(N),g=$(g)",ylims=ylims_spec,lc=:orange,lw=2,label="real type 3",guidefontsize=16,dpi=300)
plot!(ws,realparts_coex[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,label="real type 4")
plot!(ws,imagparts_coex[:,[1,2]],ylims=ylims_spec,lc=:orange,lw=2,ls=:dash,label="imag type 3")
plot!(ws,imagparts_coex[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,ls=:dash,label="imag type 4")


#plot M=2 fixed point spectrum
plot!(ws,realparts_M2[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,label="real type 3 M=2")
plot!(ws,imagparts_M2[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,ls=:dash,label="imag type 3 M=2")
plot!(ws,realparts_M2[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 M=2")
plot!(ws,imagparts_M2[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="imag type 4 M=2")

#plot M=3 fixed point spectrum
plot!(ws,realparts_M3[:,[3,4]],ylims=ylims_spec,lc=:purple,lw=2,label="real type 3 M=3")
plot!(ws,imagparts_M3[:,[3,4]],ylims=ylims_spec,lc=:purple,lw=2,ls=:dash,label="imag type 3 M=3")
plot!(ws,realparts_M3[:,[5,6]],ylims=ylims_spec,lc=:brown,lw=2,label="real type 4 M=3")
plot!(ws,imagparts_M3[:,[5,6]],ylims=ylims_spec,lc=:brown,lw=2,ls=:dash,label="imag type 4 M=3")

#plot M=4 fixed point spectrum
plot!(ws,realparts_M4[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,label="real type 3 M=4")
plot!(ws,imagparts_M4[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,ls=:dash,label="imag type 3 M=4")
plot!(ws,realparts_M4[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 M=4")
plot!(ws,imagparts_M4[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="imag type 4 M=4")

vline!([wc1],lc=:gray10,lw=2,ls=:dash,label=L"w^{(1)}_c")
vline!([wc2],lc=:gray10,lw=2,ls=:dash,label=L"w^{(2)}_c")


plot!(ws,dist_extremas_coex_mins,lc=:blue)
plot!(ws,dist_extremas_coex_maxs,lc=:blue,ls=:dash)

plot!(ws,dist_extremas_dom_mins,lc=:orange)
plot!(ws,dist_extremas_dom_maxs,lc=:orange,ls=:dash)