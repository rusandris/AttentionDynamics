using Plots
using LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("../../plotting/plot_functions/orbit_diagrams.jl")
include("../../plotting/defaults.jl")
#-----------------------plot od data-----------------------
#include standard parameters and init conditions
include("bif_floq_od_method_parameters.jl")

#create dir for figs
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

xlims=(0.5,0.65)

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
pl_od=plot(plots_od...,layout=(1,5),size=(1200,300),left_margin=8Plots.mm,bottom_margin=8Plots.mm)

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
    contres_sol_file = contres_files_N[findfirst(f -> occursin("sols",f),contres_files_N)]
    contres_bifp_file = contres_files_N[findfirst(f -> occursin("bifp",f),contres_files_N)]
   
    #read the matched files
    contres_sols = readdlm(cont_dir*contres_sol_file)
    contres_bifp = readdlm(cont_dir*contres_bifp_file)
   
    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #--------------------plot bifdiagram--------------------
    params = contres_sols[:,1]
    maxs = contres_sols[:,3]
    stability = Bool.(contres_sols[:,end])

    linestyle = [stb ? :solid : :dot for stb in stability]
    linewidth = [stb ? 2 : 1 for stb in stability]
    pl = plot(params, maxs, lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        title=L"N=%$N"*"\n"*L"g/(N-1)=%$(round(gval_phase_diag;digits=2))",
        ls=linestyle,lw=linewidth,xlabel=L"w",xlims=xlims)

    #plot some (not all) bif points

    for bp in eachrow(contres_bifp)
        type,crit_param,idx = bp
        if type == "ns"  
            scatter!([crit_param], [maxs[idx]],ms=3,mc=:orange,markerstrokewidth=0.0)
        elseif type == "nd"
            scatter!([crit_param], [maxs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
        end
    end

    #plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

    vline!(pl,[wc1],ls=:dash,lc=:gray10)
    vline!(pl,[wc2],ls=:dash,lc=:gray10)
    i == 1 && ylabel!(pl,L"max(x_1)")

    push!(branch_plots,pl)
end
pl_bif = plot(branch_plots...,layout=(1,length(Ns)),
    size=(1500,300),left_margin=8Plots.mm,bottom_margin=8Plots.mm,top_margin=8Plots.mm)

#------------------------plot Floquet exponents/multipliers--------------------------
#plot all non-unity floquet multis (real and imag parts)    
#plot trace and det? 
#plot real part?
#plot complex plane with unit circle?
plots_maxfloq = []

for (i,N) in enumerate(Ns)
    #find files that match N 
    contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]
    contres_eig_file = contres_files_N[findfirst(f -> occursin("eig",f),contres_files_N)]
    contres_sol_file = contres_files_N[findfirst(f -> occursin("sols",f),contres_files_N)]

    #find files that match N 
    od_lyap_files_N = od_lyap_files[findall(f -> occursin("N$N",f),od_lyap_files)]
    lyaps_file = od_lyap_files_N[findfirst(f -> occursin("lyaps",f),od_lyap_files_N)]

    #read the matched files
    contres_eig = readdlm(cont_dir*contres_eig_file)
    realpart = contres_eig[:,1]
    imagpart = contres_eig[:,2]
    @show maximum(abs.(imagpart))
    contres_sols = readdlm(cont_dir*contres_sol_file)
    lyaps = readdlm(od_dir*lyaps_file)[:,2]
    lyaps_params = readdlm(od_dir*lyaps_file)[:,1]

    params = contres_sols[:,1]

    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #floq maxreal
    pl_maxfloq = plot(params,realpart,mc=:gray10,xlabel=L"w",guidefontsize=guide_fontsize,
        ylims=(-1.0,0.5),xlims=xlims,markershape=:cross,ms=1,lc=:gray10,markerstrokewidth=0.0,legend=false)

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
    size=(1200,300),left_margin=8Plots.mm,bottom_margin=8Plots.mm)


pl = plot(pl_bif,pl_floq,pl_od,layout=(3,1),size=(1500,800),dpi=300)
#savefig(pl,"../../continuation_test_maxreal.png")
savefig(pl,figs_dir*"continuation_results_bif_floq_lyap_od_Ttr10000_noreinit.png")

