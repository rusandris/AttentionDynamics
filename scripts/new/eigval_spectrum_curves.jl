using Revise
#using DynamicalSystems
using AttentionDynamics
#using OrdinaryDiffEq
#using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)

#IDEA: demonstrate fixed point calclulations with figure
#make a plot where M=0, M=1, 1 < M < N and M=N cases are shown 
#plot analytical and numerical eigenvalue spectra for each case 
#discuss also g < N, and g > N cases (separate plots)

function extract_nondegenerate(vals;atol=1e-2)

    eigs_nondeg = [vals[1]]

    for i in 2:length(vals)
        val = vals[i]
        deg = any([isapprox(val, eig,atol=atol) for eig in eigs_nondeg])
        
        if !(deg)
            push!(eigs_nondeg,val)
        end
    end
    return eigs_nondeg
end

#------------------------set parameters----------------------------

#N = 4

#N = 5
#g = 5.0
N = 10 
g = 45.0 

wc1 = (N-1)/g
wc2 = (N-1)/N

#=
w_test = wc1+ (wc2-wc1)/2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
ds = create_ds(rand(2N),[g,w_test,N];diffeq=diffeq)
=#

ws = [0.0:0.01:1.0;]
#-------------------------------M = 0---------------------------------

#eigval spectrum
#coexistence fixed point 
u = zeros(2N)
numerical_spectrums = []
analytical_spectrums = []

for (i,w) in enumerate(ws)
    p = [g,w,N]

    #build Jacobian
    J = build_jacobian(u,p)

    #numerical eigenvalues
    numerical_spectrum = eigvals(J)
    vals_real = real.(numerical_spectrum)
    vals_imag = imag.(numerical_spectrum)
    vals_real_nondeg = extract_nondegenerate(vals_real)
    vals_imag_nondeg = extract_nondegenerate(vals_imag)

    push!(numerical_spectrums,vals_real_nondeg)    

    #analytical eigenvalues
    analytical_spectrum = fixed_point_eigen_spectrum(N,0,p)
    vals_real = real.(analytical_spectrum)
    vals_imag = imag.(analytical_spectrum)
    vals_real_nondeg = extract_nondegenerate(vals_real)
    vals_imag_nondeg = extract_nondegenerate(vals_imag)

    push!(analytical_spectrums,vals_real_nondeg)

end

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

#M=3 fixed point: check one of those
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=3
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
u[[2,N+2]] .= [s,s] #set both x2 and y2 to s
u[[3,N+3]] .= [s,s] #set both x3 and y3 to s

J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,3,p)

#parameter dependence -> look for bifurcations of limit cycles in terms of w

ws = [0.3:0.01:0.9;]

#reals,imags=eigen_spectrum(J,ws)
λ_specs_coex = []
λ_specs_dom = []
λ_specs_M2 = []
λ_specs_M3 = []
λ_specs_M4 = []
lyaps = []
max_crs = [] #contraction rate maxima

for (i,w) in enumerate(ws)
    p = [g,w,N]
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
    
    
    if w > wc1
        set_parameters!(ds,p)
        push!(lyaps,lyapunov(ds,1000;Ttr=100))
        traj,ts = trajectory(ds,200;Ttr=10000,Δt=1e-2)
        crs = [contraction_rate(traj[i],p) for i in 1:length(traj)]
        min_cr,max_cr = extrema(crs)
        @show max_cr
        push!(max_crs,max_cr)
    end
    
    @show w
end

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

#plot dominance fixed point spectrum 
#plot!(ws,realparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 dom")
#plot!(ws,imagparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="real type 4 dom")

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
#vline!([0.55],lc=:gray10,lw=2,ls=:dash,label=L"0.55")
zeropoint_M2 = findlast(x -> x != 0.0, imagparts_M2[:,3])
zeropoint_coex = findlast(x -> x != 0.0, imagparts_coex[:,1])

#ws[zeropoint_M2]
#ws[zeropoint_coex]
#vline!([ws[zeropoint_M2]],lc=:green,lw=2,ls=:dash,label=L"pd")
#vline!([ws[zeropoint_coex]],lc=:green,lw=2,ls=:dash,label=L"pd")

plot!(ws[25:end],lyaps,lc=:gray10,lw=2,label=L"\lambda")
plot!(ws[25:end],max_crs,lc=:purple,lw=2,label="cr max")

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(ws[end]-ws[1],6,ws[1],-6), opacity=.2,fillc=:gray50,label="")

savefig(pl,"../../eigen_spectrum_coexistence_fp_N$(N).png")
