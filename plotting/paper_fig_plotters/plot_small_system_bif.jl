using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
using LinearAlgebra
cd(@__DIR__)
include("defaults.jl")

figs_dir = "../../paper_figs/"
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
cont_dir = "../../data/bifkit_data_g5_N2/"
#cont_dir = "../../bifkit_data_gN5/"


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
contres_sol_file = contres_files_N[findfirst(f -> occursin("contres_po_sols",f),contres_files_N)]
contres_bifp_file = contres_files_N[findfirst(f -> occursin("contres_po_bifp",f),contres_files_N)]

#read the matched files
contres_sols = readdlm(cont_dir*contres_sol_file)
contres_bifp = readdlm(cont_dir*contres_bifp_file)


#--------------------plot bifdiagram------------------
lc_color = collect(Colors.JULIA_LOGO_COLORS)[[1]]  #use logo colors

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

annotate!(pl_bif,0.15, 4.15, text(L"\textnormal{coexistence}", :gray10, :center, 15))
annotate!(pl_bif,0.55, 4.15, text(L"\textnormal{dominance}", :gray10, :center, 15))

#--------------------------plot critical values------------------------

rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
plot!(rectangle(0.1, 4, 0.1, 0), opacity=0.2, c=:gray10)
plot!(rectangle(0.1, 4, 0.5, 0), opacity=0.2, c=:gray10)


vline!(pl_bif, [wc1], ls=:dash, lw=3, lc=:gray15)
vline!(pl_bif, [wc2], ls=:dot, lw=3, lc=:gray15)

#w_c^{(1)} = \epsilon (N-1)
#w_c^{(2)} = \frac{N-1}{N}
annotate!(wc1, 0.0, text(L"w_{\mathrm{c}}", :gray10, :center, 22))
annotate!(wc2, 0.0, text(L"w_{\mathrm{d}}", :gray10, :center, 22))

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

params = contres_sols[:,1]
maxs = contres_sols[:,3]
mins = contres_sols[:,2]
stability = Bool.(contres_sols[:,end])

linestyle = [stb ? :solid : :dot for stb in stability]
linewidth = [stb ? 4 : 2 for stb in stability]
#sol_sampling = [stb ? 1 : 20 for stb in stability]

plot!(pl_bif,params[1:plot_srate:end], maxs[1:plot_srate:end],lc=lc_color,legend=false,
ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end])

#annotate!(pl_bif,0.35, 2.5, text(L"\textnormal{stable limit cycle}", :gray10, :center, 15))

#=

#plot stable and unstable parts separately
maxs_stable = maxs[stability]  #select stable 
maxs_unstable = maxs[.!(stability)] #negation

plot!(pl,params[stability][1:plot_srate:end], maxs_stable[1:plot_srate:end], lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=:solid,lw=1,xformatter=:auto,xlims=xlims_bif_pl)

plot!(pl,params[.!(stability)][1:plot_srate:end], maxs_unstable[1:plot_srate:end], mc=:gray10,legend=false,guidefontsize=guide_fontsize,
    st=:scatter,ms=1,xformatter=:auto,xlims=xlims_bif_pl)

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


#-------------------------------bif points----------------------------
annotate!(0.22, 1.15, text(L"H", :gray10, :center, 18))
annotate!(0.48, 2.2, text(L"HC", :gray10, :center, 18))

annotate!(0.15,0.7, text(L"\mathbf{u}^{*}_{\mathrm{c}}", :gray10, :center, 20))
annotate!(0.55,1.9, text(L"\mathbf{u}^{*}_{\mathrm{d}}", :gray10, :center, 20))

#--------------------------plot inset----------------------------
colors = collect(Colors.JULIA_LOGO_COLORS)[[1, 2]]  #use logo colors

#--------------------coexistence---------------------
set_parameter!(ds,2,0.15)
traj1,t = trajectory(ds,10;Ttr=0)
plot!(inset = (1, bbox(0.06,0.1,0.3,0.2)), subplot = 2)
plot!(pl_bif[2],t,traj1[:,1],lw=3,
    label=L"x_1",lc=colors[1],
    yticks=[0,2,4],ylims=(0.0,4), xlims=(-0.5,10),
    xticks=[0,5,10],tickfontsize=tick_fontsize_small,
    legend=false)
    #legend_columns=-1, legend=:topright,
    #legendfontsize=legend_fontsize-7)
plot!(pl_bif[2], t, traj1[:, 3], lw=3, label="", lc=colors[2])
annotate!(pl_bif[2],-2, 2.0, text(L"x_i", :gray10, :center, 14))
annotate!(pl_bif[2], 5, -1.6, text(L"t", :gray10, :center, 14))
plot!(pl_bif[1], [0.15, 0.18], [1.1, 2.8], arrow=true, lcolor=:gray10, linewidth=2, label="")


hline!(pl_bif[2], [1.0], lw=2, lc=:gray10, ls=:dash, xlabel="", label="")


#--------------------WLC---------------------
set_parameter!(ds, 2, 0.35)
traj0, t = trajectory(ds, 40; Ttr=100)
plot!(inset=(1, bbox(0.35, 0.5, 0.3, 0.2)), subplot=3)
plot!(pl_bif[3], t, traj0[:, 1], lw=3,
    label=L"x_1", lc=colors[1],
    yticks=[0, 2, 4], ylims=(0, 4),
    xticks=[100, 120, 140], tickfontsize=tick_fontsize_small,
    legend=false)
#legend_columns=-1, legend=:topright,
#legendfontsize=legend_fontsize-7)
plot!(pl_bif[3], t, traj0[:, 3], lw=3, label="", lc=colors[2])
annotate!(pl_bif[3], 92, 2.0, text(L"x_i", :gray10, :center, 14))
annotate!(pl_bif[3], 120, -1.6, text(L"t", :gray10, :center, 14))
plot!(pl_bif[1],[0.35, 0.35], [3, 2.5], arrow=true, color=lc_color, linewidth=2, label="")

hline!(pl_bif[2], [1.0], lw=2, lc=:gray10, ls=:dash, xlabel="", label="")


#--------------------dominance---------------------
set_parameter!(ds,2,0.55)
traj2,t = trajectory(ds,10;Ttr=0)
plot!(inset = (1, bbox(0.70,0.1,0.3,0.2)), subplot = 4)
plot!(pl_bif[4],t,traj2[:,1],lw=3,
    label="", lc=colors[1],
    yticks=[0,2,4], xlims=(-0.5, 10), ylims=(0, 4),
    xticks=[0, 5, 10], tickfontsize=tick_fontsize_small,
    legend=false)
    #legend=(0.45, 0.28),legend_columns=-1,legendfontsize=legend_fontsize-7)
plot!(pl_bif[4], t, traj2[:, 3], lw=3, label=L"x_2", lc=colors[2])
annotate!(pl_bif[4], -2, 2.0, text(L"x_i", :gray10, :center, 14))
annotate!(pl_bif[4], 5, -1.6, text(L"t", :gray10, :center, 14))
hline!(pl_bif[4], [1 / (1 - 0.55)], lw=2, lc=:gray10, ls=:dash, label="") #xlabel=L"t")
plot!(pl_bif[1], [0.55, 0.55], [2.5, 2.8], arrow=true, color=:gray10, linewidth=2, label="")

#add title 
plot!(pl_bif,plot_title=L"N=2",size=(700,400))

#add plot label annotation 
annotate!(pl_bif[1], 0.07, 4.3, text(L"(a)", :gray10, :center, 18))

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


