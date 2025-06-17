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

xlims = (ws[1],ws[end])
ylims_bif = (0.0,4.0)

#plot coexistence fixed point spectrum
pl_eig = plot(ws,realparts_coex[:,[1,2]],ylabel="",xlabel=L"w",
    legend=false,
    #markershape=:circle,
    #title="N=$(N),g=$(g)",
    mc=:orange,ms=3,
    xlims=xlims,ylims=[-1.0,4.0],
    lc=:orange,
    label="",guidefontsize=16,dpi=300)
plot!(pl_eig,ws,imagparts_coex[:,[1,2]],lc=:orange,ms=3,lw=2,label="")
plot!(pl_eig,ws,realparts_dom[:,[2]],lc=:blue,ms=3,mc=:blue,lw=2,label="")
hline!(pl_eig,[0.0],ls=:dash,lw=2,lc=:gray10,)

#plot critical values 
vline!(pl_eig,[wc1],ls=:dash,lw=2,lc=:gray10)
vline!(pl_eig,[wc2],ls=:dash,lw=2,lc=:gray10)

#plot dominance fixed point spectrum 
#plot!(ws,realparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 dom")
#plot!(ws,imagparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="real type 4 dom")

#--------------------lyapunov--------------------------
d0 = 1e-9
states = [u0,u0 .+ d0/sqrt(length(u0))]
pds = ParallelDynamicalSystem(ds,states)
@show N
lyaps = zeros(length(ws))
for (i,w) in enumerate(ws)
    @show w
    @show current_parameters(pds)
    #@show current_states(pds)
    reinit!(pds)
    set_parameter!(pds,2,w)
    lyaps[i] = lyapunov(pds,100;Ttr=200,Δt=1.0)
    @show lyaps[i]
end
plot!(pl_eig,ws,lyaps,lw=5,lc=:gray10)


#------------------------bifurcation diagram---------------------------------

#data dirs
cont_dir = "../../plotting/DDE25/bifkit_data_g5_N2/"
#cont_dir = "../../bifkit_data_gN5/"

#create dir for figs
figs_dir = "./poster_figs/"
mkpath(figs_dir)
plot_srate = 1

#saved params
global_params = readdlm(cont_dir*"global_params_Ns_g.txt")
Ns = Int.(global_params[1:end-1])
N = Ns[1]
g = global_params[end]

contres_files = readdir(cont_dir)


#read cont results

#find files that match N 
contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]

#limit cycle continuation results
contres_sol_file = contres_files_N[findfirst(f -> occursin("contres_po_sols",f),contres_files_N)]
contres_bifp_file = contres_files_N[findfirst(f -> occursin("contres_po_bifp",f),contres_files_N)]

#read the matched files
contres_sols = readdlm(cont_dir*contres_sol_file)
contres_bifp = readdlm(cont_dir*contres_bifp_file)

#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N

#--------------------plot bifdiagram------------------

#-------------------plot fixed point results----------------------
#use analytical results for fixed points
xs_fp = ones(length(ws))

stability_fp = Bool.( ws .< wc1 )
stability_fp2 = Bool.( ws .> wc2 )

linestyle_fp = [stb ? :solid : :dash for stb in stability_fp]
linewidth_fp = [stb ? 2 : 1 for stb in stability_fp]

linestyle_fp2 = [stb ? :solid : :dash for stb in stability_fp2]
linewidth_fp2 = [stb ? 2 : 1 for stb in stability_fp2]


pl_bif = plot(ws, xs_fp, lc=:gray10,
    legend=false,guidefontsize=guide_fontsize,xlims=xlims,ylims=ylims_bif,
    #title=L"N=%$N"*"\n"*L"g/(N-1)=%$(round(gval_phase_diag;digits=2))",
    ls=linestyle_fp,lw=linewidth_fp,xformatter=:none)

fp2(w) = 1/(1-w)
xs_fp2 = fp2.(ws)

plot!(ws, xs_fp2, lc=:gray10,
    legend=false,guidefontsize=guide_fontsize,xlims=xlims,ylims=ylims_bif,
    ls=linestyle_fp2,lw=linewidth_fp2)

annotate!(pl_bif,0.15, 3.0, text("stable fixed point \n (coexistence)", :gray10, :center, 6))
annotate!(pl_bif,0.55, 3.0, text("stable fixed point \n (dominance)", :gray10, :center, 6))
scatter!([wc1],[1.0],ms=3,mc=:gray10) #coexistence
scatter!([wc2],[fp2(wc2)],ms=3,mc=:gray10) #dominance

#=

#plot some (not all) bif points
for bp in eachrow(contres_fp_bifp)
    type,crit_param,idx = bp
    if type == "ns"  
        scatter!([crit_param], [xs[idx]],ms=3,mc=:orange,markerstrokewidth=0.0)
    elseif type == "nd"
        scatter!([crit_param], [xs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
    end
end
=#

#-------------------plot lim cycle results----------------------

params = contres_sols[:,1]
maxs = contres_sols[:,3]
mins = contres_sols[:,2]
stability = Bool.(contres_sols[:,end])

linestyle = [stb ? :solid : :dot for stb in stability]
linewidth = [stb ? 2 : 1 for stb in stability]
#sol_sampling = [stb ? 1 : 20 for stb in stability]

plot!(pl_bif,params[1:plot_srate:end], maxs[1:plot_srate:end],lc=:gray10,legend=false,guidefontsize=guide_fontsize,
ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end])

annotate!(pl_bif,0.35, 3.5, text("stable limit cycle", :gray10, :center, 8))

#=

#plot stable and unstable parts separately
maxs_stable = maxs[stability]  #select stable 
maxs_unstable = maxs[.!(stability)] #negation

plot!(pl,params[stability][1:plot_srate:end], maxs_stable[1:plot_srate:end], lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=:solid,lw=1,xformatter=:auto,xlims=xlims_pl)

plot!(pl,params[.!(stability)][1:plot_srate:end], maxs_unstable[1:plot_srate:end], mc=:gray10,legend=false,guidefontsize=guide_fontsize,
    st=:scatter,ms=1,xformatter=:auto,xlims=xlims_pl)

=#

#plot some (not all) bif points

for bp in eachrow(contres_bifp)
    type,crit_param,idx = bp
    if type == "ns"  
        scatter!([crit_param], [maxs[idx]],ms=3,mc=:orange,markerstrokewidth=0.0)
    elseif type == "nd"
        scatter!([crit_param], [maxs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
    elseif type == "bp"
        scatter!([crit_param], [maxs[idx]],ms=3,mc=:green,markerstrokewidth=0.0)
    end
end
#plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

vline!(pl_bif,[wc1],ls=:dash,lw=2,lc=:gray10)
vline!(pl_bif,[wc2],ls=:dash,lw=2,lc=:gray10)


l = @layout [a{0.8h}; b{0.2h}]
pl = plot(pl_bif,pl_eig,layout=l,size=(500,300))

savefig(pl,figs_dir*"system_small_N2_bifdiag.png")
