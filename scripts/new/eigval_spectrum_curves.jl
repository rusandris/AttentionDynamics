using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)

#set parameters
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
g_N = gs[30]

#N = 5
#N = 4
N = 7 
g = g_N * (N-1) 
wc1 = (N-1)/g
wc2 = (N-1)/N
w_test = wc1+ (wc2-wc1)/2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)

ds = create_ds(rand(2N),[g,w_test,N];diffeq=diffeq)

#eigval spectrum
#coexistence fixed point 
u = ones(2N)
p = current_parameters(ds)
J = build_jacobian(u,p)

numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,N,p)
#ws = [wc1+0.02:0.001:wc2-0.01;] 
ws = [0.3:0.001:0.9;]

#reals,imags=eigen_spectrum(J,ws)
λ_specs = []
for (i,w) in enumerate(ws)
    p = [g,w,N]
    spec = fixed_point_eigen_spectrum(N,N,p)
    push!(λ_specs,spec)
end

realparts = stack(real.(λ_specs))'
imagparts = stack(imag.(λ_specs))'

ylims_spec = [-6.0,3.0]
pl = plot(ws,realparts[:,[1,2]],ylabel="eigen spectrum",xlabel=L"w",
    title="N=$(N),g=$(g)",ylims=ylims_spec,lc=:orange,lw=2,label="real type 3",guidefontsize=16,dpi=300)
plot!(ws,realparts[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,label="real type 4")
plot!(ws,imagparts[:,[1,2]],ylims=ylims_spec,lc=:orange,lw=2,ls=:dash,label="imag type 3")
plot!(ws,imagparts[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,ls=:dash,label="imag type 4")

vline!([wc1],lc=:gray10,lw=2,ls=:dash,label=L"w^{(1)}_c")
vline!([wc2],lc=:gray10,lw=2,ls=:dash,label=L"w^{(2)}_c")
#vline!([0.55],lc=:gray10,lw=2,ls=:dash,label=L"0.55")

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(ws[end]-ws[1],6,ws[1],-6), opacity=.2,fillc=:gray50,label="")

savefig(pl,"../../eigen_spectrum_coexistence_fp_N$(N).png")
