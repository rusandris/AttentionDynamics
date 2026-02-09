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


N = 2
g = 5
wc2 = (N-1)/N 
wc1 = (N-1)/g 
w = 0.1

p = [g,w,N]


ws = [0.1:0.005:0.6;]
#critical values
wc1 = (N-1)/g
wc2 = (N-1)/N

#u0 = generate_initial_condition(N)
u0 = [0.5,0.0,1.3,0.0] 
#u0 = rand(2N)
ds = create_ds(u0,p)


#------------------------bifurcation diagram---------------------------------

#data dirs
cont_dir = "../../data/bifkit_data_g5_N2_trap/"
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

xlims_bif = (ws[1],ws[end])
ylims_bif = (0.5,4.1)


#read cont results

#find files that match N 
contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]

#limit cycle continuation results
contres_sol_files = contres_files_N[findall(f -> occursin("contres_po_sols",f),contres_files_N)]
contres_bifp_files = contres_files_N[findall(f -> occursin("contres_po_bifp",f),contres_files_N)]




#--------------------plot bifdiagram------------------

#-------------------plot fixed point results----------------------
#use analytical results for fixed points
xs_fp = ones(length(ws))

stability_fp = Bool.( ws .<= wc1 )
stability_fp2 = Bool.( ws .>= wc2 )

linestyle_fp = [stb ? :solid : :dot for stb in stability_fp]
linewidth_fp = [stb ? 4 : 2 for stb in stability_fp]

linestyle_fp2 = [stb ? :solid : :dot for stb in stability_fp2]
linewidth_fp2 = [stb ? 4 : 2 for stb in stability_fp2]


pl_bif = plot(ws, xs_fp, lc=:gray10,
    legend=false,
    xlims=xlims_bif,ylims=ylims_bif,
    tickfontsize=tick_fontsize_large,
    guidefontsize=guide_fontsize_large,
    xlabel=L"w",
    ylabel=L"\mathrm{max}(x_1)",
    #title=L"N=%$N"*"\n"*L"g/(N-1)=%$(round(gval_phase_diag;digits=2))",
    ls=linestyle_fp,lw=linewidth_fp,
    size=(800,500),
    xformatter=:auto,
    right_margin=6Plots.mm,
    dpi=300)

fp2(w) = 1/(1-w)
xs_fp2 = fp2.(ws)

plot!(ws, xs_fp2, lc=:gray10,
    legend=false,xlims=xlims_bif,ylims=ylims_bif,
    ls=linestyle_fp2,lw=linewidth_fp2)

annotate!(pl_bif,0.15, 3.8, text("coexistence", :gray10, :center, 15))
annotate!(pl_bif,0.55, 3.8, text("dominance", :gray10, :center, 15))


#-------------------------plot bif points--------------------------
scatter!([wc1],[1.0],ms=9,mc=:gray10) #coexistence
scatter!([wc2],[fp2(wc2)],ms=9,mc=:gray10) #dominance


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
for i in 1:length(contres_sol_files)

    contres_sol_file = contres_sol_files[i]
    contres_bifp_file = contres_sol_files[i]

    #read the matched files
    contres_sols = readdlm(cont_dir * contres_sol_file)
    contres_bifp = readdlm(cont_dir * contres_bifp_file)

    params = contres_sols[:,1]
    maxs = contres_sols[:,3]
    mins = contres_sols[:,2]
    stability = Bool.(contres_sols[:,end])

    linestyle = [stb ? :solid : :dot for stb in stability]
    linewidth = [stb ? 4 : 2 for stb in stability]
    #sol_sampling = [stb ? 1 : 20 for stb in stability]

    plot!(pl_bif,params[1:plot_srate:end], maxs[1:plot_srate:end],lc=:gray10,legend=false,
    ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end])

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
    annotate!(pl_bif,0.35, 2.5, text("stable limit cycle", :gray10, :center, 15))
end

#=

#plot stable and unstable parts separately
maxs_stable = maxs[stability]  #select stable 
maxs_unstable = maxs[.!(stability)] #negation

plot!(pl,params[stability][1:plot_srate:end], maxs_stable[1:plot_srate:end], lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=:solid,lw=1,xformatter=:auto,xlims=xlims_bif_pl)

plot!(pl,params[.!(stability)][1:plot_srate:end], maxs_unstable[1:plot_srate:end], mc=:gray10,legend=false,guidefontsize=guide_fontsize,
    st=:scatter,ms=1,xformatter=:auto,xlims=xlims_bif_pl)

=#


#plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

#--------------------------plot critical values------------------------

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(0.1,4,0.1,0), opacity=.2,c=:gray10)
plot!(rectangle(0.1,4,0.5,0), opacity=.2,c=:gray10)

# add vertical lines 

vline!(pl_bif,[wc1],ls=:dash,lw=3,lc=:gray10)
vline!(pl_bif,[wc2],ls=:dash,lw=3,lc=:red)
vline!(pl_bif, [0.3], ls=:dash, lw=1, lc=:gray50)
vline!(pl_bif, [0.41], ls=:dash, lw=1, lc=:gray50)


#w_c^{(1)} = \epsilon (N-1)
#w_c^{(2)} = \frac{N-1}{N}
annotate!(wc1+0.03, 3.8, text(L"w_c^{(1)}", :gray10, :center, 22))
annotate!(wc2-0.03, 3.8, text(L"w_c^{(2)}", :gray10, :center, 22))


#-------------------------------bif points----------------------------
annotate!(0.22, 0.8, text("HOPF", :gray10, :center, 13))
annotate!(0.45, 1.5, text("HOMOCLINIC", :gray10, :center, 13))

annotate!(0.15,0.8, text(L"u^{*}_N", :gray10, :center, 22))
annotate!(0.55,1.8, text(L"u_1^{*}", :gray10, :center, 22))

#--------------------------plot inset----------------------------

#--------------------coexistence---------------------
set_parameter!(ds,2,0.15)
traj1,t = trajectory(ds,10;Ttr=0)
plot!(inset = (1, bbox(0.05,0.15,0.25,0.25)), subplot = 2)
hline!(pl_bif[2],[1.0],lw=2,lc=:gray10,ls=:dash,xlabel="",label="")
plot!(pl_bif[2],t,traj1[:,1],lw=3,label=L"x_1",
    legend_columns=-1,legend=:bottom,legendfontsize=legend_fontsize-5)
plot!(pl_bif[2],t,traj1[:,2],lw=3,label=L"x_2")
annotate!(pl_bif[2],10, 0.15, text(L"t", :gray10, :center, 12))


#--------------------dominance---------------------
set_parameter!(ds,2,0.55)
traj2,t = trajectory(ds,10;Ttr=0)
plot!(inset = (1, bbox(0.74,0.15,0.25,0.25)), subplot = 3)
hline!(pl_bif[3],[1/(1-0.55)],lw=2,lc=:gray10,ls=:dash,label="") #xlabel=L"t")
plot!(pl_bif[3],t,traj2[:,1],lw=3,legend=(0.45,0.28),label=L"x_1",
    legend_columns=-1,legendfontsize=legend_fontsize-5)
plot!(pl_bif[3],t,traj2[:,3],lw=3,label=L"x_2")
annotate!(pl_bif[3],10, 0.3, text(L"t", :gray10, :center, 12))


pl_bif
savefig(pl_bif,figs_dir*"system_small_N2_bifdiag_complicated.pdf")


#l = @layout [a{0.8h}; b{0.2h}]
#pl = plot(pl_bif,pl_eig,layout=l,size=(500,300))

#savefig(pl,figs_dir*"system_small_N2_bifdiag.png")


#---------------------------------eigvals,lyap----------------------------------------
#=

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
=#

#=
#---------------------calculate eigenvalues for a range of ws--------------------


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
pl_eig = plot(ws,realparts_coex[:,[1,2]],ylabel="",xlabel=L"w",
    legend=false,
    #markershape=:circle,
    #title="N=$(N),g=$(g)",
    mc=:orange,ms=3,
    xlims=xlims_bif,ylims=[-1.0,4.0],
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
=#


