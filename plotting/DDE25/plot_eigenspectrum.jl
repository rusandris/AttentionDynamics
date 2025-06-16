using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
using LinearAlgebra
cd(@__DIR__)
include("../defaults.jl")

figs_dir = "./poster_figs/"
mkpath(figs_dir)

logocolors = Colors.JULIA_LOGO_COLORS
ms=6

function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end

N = 2
g = 5
wc2 = (N-1)/N 
wc1 = (N-1)/g 
w = 0.1

p = [g,w,N]

#u0 = generate_initial_condition(N)
u0 = [0.5,0.0,1.3,0.0] 
#u0 = rand(2N)
ds = create_ds(u0,p)

#eigval spectrum
#coexistence fixed point 
u = ones(2N)
p = current_parameters(ds)
J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,N,p)

#dominance fixed point: doesn't matter which one we choose 
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=1
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)

u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,1,p)

#M=2 fixed point: check one of those
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=2
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
u[[2,N+2]] .= [s,s] #set both x2 and y2 to s

J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,2,p)


#---------------------calculate eigenvalues for a range of ws--------------------
ws = [0.1:0.01:0.6;]

#reals,imags=eigen_spectrum(J,ws)
λ_specs_coex = []
λ_specs_dom = []

for (i,w) in enumerate(ws)
    p = [g,w,N]
    spec_coex = fixed_point_eigen_spectrum(N,N,p) #coexistence spectrum
    spec_dom = fixed_point_eigen_spectrum(N,1,p) #dominance spectrum

    push!(λ_specs_coex,spec_coex)
    push!(λ_specs_dom,spec_dom)
end

#extract real and imaginary parts 
realparts_coex = stack(real.(λ_specs_coex))'
imagparts_coex = stack(imag.(λ_specs_coex))'
realparts_dom = stack(real.(λ_specs_dom))'
imagparts_dom = stack(imag.(λ_specs_dom))'

#plot coexistence fixed point spectrum
pl = plot(ws,realparts_coex[:,[1,2]],ylabel="eigen spectrum",xlabel=L"w",
    legend=false,
    title="N=$(N),g=$(g)",lc=:orange,lw=2,
    size=(800,300),
    label="",guidefontsize=16,dpi=300)
plot!(ws,realparts_dom[:,[2]],lc=:blue,lw=2,label="")
hline!([0.0],ls=:dash,lw=2,lc=:gray10,)

#plot critical values 
vline!([wc1],ls=:dash,lw=2,lc=:gray10)
vline!([wc2],ls=:dash,lw=2,lc=:gray10)

#plot dominance fixed point spectrum 
#plot!(ws,realparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 dom")
#plot!(ws,imagparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="real type 4 dom")