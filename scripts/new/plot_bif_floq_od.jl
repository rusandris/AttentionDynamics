using Plots
using LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("../../plotting/plot_functions/orbit_diagrams.jl")
include("../../plotting/defaults.jl")
#-----------------------plot od data-----------------------
#include standard parameters and init conditions
include("bif_floq_od_method_parameters.jl")

#data dirs
od_dir = "../../orbit_diagram_data_gN5/"
cont_dir = "../../bifkit_data_gN5_210/"
#cont_dir = "../../bifkit_data_gN5/"

#create dir for figs
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

#xlims=(0.5,0.65)
xlims=(0.1,0.4)
plot_srate = 1

od_lyap_files = readdir(od_dir)
plots_od = []
for (i,N) in enumerate(Ns)

    od_lyap_files_N = od_lyap_files[findall(f -> occursin("N$N",f),od_lyap_files)]
    od_file = od_lyap_files_N[findfirst(f -> occursin("od",f),od_lyap_files_N)]

    #read the matched files
    od_data = readdlm(od_dir*od_file)

    ps = od_data[1,:]   
    od = od_data[2:end,:]

    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    pl = plot_orbit_diagram(ps,od;OD_DEFAULTS...,ma=0.4,ms=0.5,tickfontsize=8,ylabel="",xlims=xlims,widen=false)
    i == 1 && ylabel!(pl,L"x_1 pmap")

    vline!([wc1],ls=:dash,lc=:gray10)
    vline!([wc2],ls=:dash,lc=:gray10)

    push!(plots_od,pl)
end
pl_od=plot(plots_od...,layout=(1,5),size=(1200,300),left_margin=8Plots.mm,bottom_margin=5Plots.mm);

#----------------plot each branch of periodic orbits---------------------

#saved params
global_params = readdlm(cont_dir*"global_params_Ns_g.txt")
Ns = Int.(global_params[1:end-1])
gval_phase_diag = global_params[end]

contres_files = readdir(cont_dir)

branch_plots = []
for (i,N) in enumerate(Ns)
    #read cont results
    
    #find files that match N 
    contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]

    #fixed point continuation results
    contres_fp_sol_file = contres_files_N[findfirst(f -> occursin("contres_fp_sols",f),contres_files_N)]
    contres_fp_bifp_file = contres_files_N[findfirst(f -> occursin("contres_fp_bifp",f),contres_files_N)]

    #read the matched files
    contres_fp_sols = readdlm(cont_dir*contres_fp_sol_file)
    contres_fp_bifp = readdlm(cont_dir*contres_fp_bifp_file)

    #limit cycle continuation results
    contres_sol_file = contres_files_N[findfirst(f -> occursin("contres_po_sols",f),contres_files_N)]
    contres_bifp_file = contres_files_N[findfirst(f -> occursin("contres_po_bifp",f),contres_files_N)]
   
    #read the matched files
    contres_sols = readdlm(cont_dir*contres_sol_file)
    contres_bifp = readdlm(cont_dir*contres_bifp_file)
   
    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #--------------------plot bifdiagram------------------

    params = contres_fp_sols[:,1]
    xs = contres_fp_sols[:,2]
    @show extrema(xs)
    stability = Bool.(contres_fp_sols[:,end])
    @show extrema(stability)

    linestyle = [stb ? :solid : :dash for stb in stability]
    linewidth = [stb ? 1 : 1 for stb in stability]

    #-------------------plot fixed point results----------------------
    pl = plot(params, xs, lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        title=L"N=%$N"*"\n"*L"g/(N-1)=%$(round(gval_phase_diag;digits=2))",
        ls=linestyle,lw=linewidth,xformatter=:none,xlims=xlims)

    display(pl)
    @show length(contres_fp_bifp)
    #plot some (not all) bif points
    for bp in eachrow(contres_fp_bifp)
        type,crit_param,idx = bp
        if type == "ns"  
            scatter!([crit_param], [xs[idx]],ms=3,mc=:orange,markerstrokewidth=0.0)
        elseif type == "nd"
            scatter!([crit_param], [xs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
        end
    end

    #-------------------plot lim cycle results----------------------

    params = contres_sols[:,1]
    maxs = contres_sols[:,3]
    mins = contres_sols[:,2]
    @show extrema(maxs)
    stability = Bool.(contres_sols[:,end])
    @show extrema(stability)

    linestyle = [stb ? :solid : :dot for stb in stability]
    linewidth = [stb ? 2 : 1 for stb in stability]
    #sol_sampling = [stb ? 1 : 20 for stb in stability]

    plot!(pl,params[1:plot_srate:end], maxs[1:plot_srate:end],lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end],xformatter=:none,xlims=xlims)

    #=

    #plot stable and unstable parts separately
    maxs_stable = maxs[stability]  #select stable 
    maxs_unstable = maxs[.!(stability)] #negation

    plot!(pl,params[stability][1:plot_srate:end], maxs_stable[1:plot_srate:end], lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        ls=:solid,lw=1,xformatter=:none,xlims=xlims)

    plot!(pl,params[.!(stability)][1:plot_srate:end], maxs_unstable[1:plot_srate:end], mc=:gray10,legend=false,guidefontsize=guide_fontsize,
        st=:scatter,ms=1,xformatter=:none,xlims=xlims)

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
    display(pl)
    #plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

    vline!(pl,[wc1],ls=:dash,lc=:gray10)
    vline!(pl,[wc2],ls=:dash,lc=:gray10)
    i == 1 && ylabel!(pl,L"max(x_1)")
    display(pl)
    push!(branch_plots,pl)
end
pl_bif = plot(branch_plots...,layout=(1,length(Ns)),
    size=(1500,300),left_margin=8Plots.mm)#bottom_margin=reduced_bottom_margin)

#------------------------plot Floquet exponents/multipliers--------------------------
#plot all non-unity floquet multis (real and imag parts)    
#plot trace and det? 
#plot real part?
#plot complex plane with unit circle?
plots_maxfloq = []

for (i,N) in enumerate(Ns)
    #find files that match N 
    contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]
    contres_eig_file = contres_files_N[findfirst(f -> occursin("contres_po_eig",f),contres_files_N)]
    contres_sol_file = contres_files_N[findfirst(f -> occursin("contres_po_sols",f),contres_files_N)]

    #find files that match N 
    od_lyap_files_N = od_lyap_files[findall(f -> occursin("N$N",f),od_lyap_files)]
    lyaps_file = od_lyap_files_N[findfirst(f -> occursin("lyaps",f),od_lyap_files_N)]

    #read the matched files
    contres_eig = readdlm(cont_dir*contres_eig_file)

    #filter eig results
    #contains real and imag parts separately
    #

    realparts = contres_eig[:,1:2N]
    imagparts = contres_eig[:,(2N+1):end]

    #how to filter out zero exponents
    
    realparts_max = maximum(realparts;dims=2)


    contres_sols = readdlm(cont_dir*contres_sol_file)
    lyaps = readdlm(od_dir*lyaps_file)[:,2]
    lyaps_params = readdlm(od_dir*lyaps_file)[:,1]

    params = contres_sols[:,1]

    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #floq maxreal
    @show length(params),length(realparts_max)
    pl_maxfloq = plot(params,realparts,mc=:gray10,xformatter=:none,guidefontsize=guide_fontsize,
        ylims=(Inf,Inf),xlims=xlims,markershape=:cross,ms=1,lc=:gray10,markerstrokewidth=0.0,legend=false)

    #floq maximag
    #plot!(pl_maxfloq,params,imagpart,mc=:gray50,xlabel=L"w",guidefontsize=guide_fontsize,
    #    xlims=xlims,markershape=:cross,ms=1,lc=:gray50,markerstrokewidth=0.0,legend=false)

    #lyapunov LLE
    plot!(lyaps_params,lyaps,mc=:gray10,markershape=:circle,ms=1,lc=:red,markerstrokewidth=0.0,legend=false)

    vline!([wc1],ls=:dash,lc=:gray10)
    vline!([wc2],ls=:dash,lc=:gray10)

    hline!([0.0],ls=:dash)
    i == 1 && ylabel!(pl_maxfloq,L"max(Re(μ_i)),\lambda")

    push!(plots_maxfloq,pl_maxfloq)

end
pl_floq = plot(plots_maxfloq...,layout=(1,length(Ns)),
    size=(1200,300),left_margin=8Plots.mm)#bottom_margin=reduced_bottom_margin)


pl = plot(pl_bif,pl_floq,pl_od,layout=(3,1),size=(1500,600),dpi=300)
#savefig(pl,"../../continuation_test_maxreal.png")
savefig(pl,figs_dir*"continuation_results_bif_floq_lyap_od_Ttr20000_noreinit_g5_N26.png")

