using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("..//defaults.jl")

dir_cont = "../../data/bifkit_data_orig_N5_around_hopf/"

#---------------plotting params---------------

#xlims=(0.5,0.65)
xlims_pl=(0.18,0.4)
plot_srate = 1
#od_lyap_files = readdir(od_dir)
plots_od = []
tickfontsize_od_bif = 10
yticks_od = [[0:0.5:1;],[0:0.5:1;],[0:5:15;]]
ylims_od = [(0,1.1),(0,1.2),(0,15)]
reduced_margin = -9Plots.mm
reduced_margin_side = -2Plots.mm
ms_bp = 5 
annotation_size = 12



#create dir for figs
figs_dir = "../paper_figs/poster_figs/"
#mkpath(figs_dir)

#--------------read continuation results-------------------
#read all attractor continuation
#plot on the same figure 
dir_attrs = readdir(dir_cont)

pl = plot(; legend=false,guidefontsize=guide_fontsize,
    xformatter=:auto,widen=true,
    DEFAULTS...,tickfontsize=tick_fontsize+3,yformatter=:auto)


for (i,dir_attr) in enumerate(dir_attrs)
    #saved params
    @show i

    dir_cont_i = dir_cont*dir_attr*"/"
    global_params = readdlm(dir_cont_i*"global_params_Ns_g.txt")
    N = Int.(global_params[1])
    g = global_params[2]
    c_crit1 = global_params[end-1]
    c_crit2 = global_params[end]

    #read cont results
    contres_files = readdir(dir_cont_i)
    branch_plots = []
    @show contres_files

    #find files that match N 
    contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]

    #limit cycle continuation results
    contres_sol_files = contres_files_N[findall(f -> occursin("contres_po_sols",f),contres_files_N)]
    contres_bifp_files = contres_files_N[findall(f -> occursin("contres_po_bifp",f),contres_files_N)]

    #read the matched files
    contres_sols = readdlm.(dir_cont_i .* contres_sol_files)
    contres_bifp = readdlm.(dir_cont_i .* contres_bifp_files)

    #-------------------plot lim cycle results----------------------

    limc_sols = contres_sols[1]
    limc_bifp = contres_bifp[1]

    params = limc_sols[:,1]
    @show extrema(params)
    maxs = limc_sols[:,3]
    mins = limc_sols[:,2]
    stability = Bool.(limc_sols[:,end])

    linestyle = [stb ? :solid : :dot for stb in stability]
    linewidth = [stb ? 4 : 1 for stb in stability]
    #sol_sampling = [stb ? 1 : 20 for stb in stability]

    plot!(pl,params[1:plot_srate:end], maxs[1:plot_srate:end],
        lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end])

    #plot some (not all) bif points

    for bp in eachrow(limc_bifp)
        type,crit_param,idx = bp
        if type == "ns"  
            scatter!([crit_param], [maxs[idx]],ms=ms_bp,mc=:orange,markerstrokewidth=0.0)
        elseif type == "nd"
            scatter!([crit_param], [maxs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
        elseif type == "bp"
            scatter!([crit_param], [maxs[idx]],ms=3,mc=:green,markerstrokewidth=0.0)
        end
    end

    vline!([c_crit1],ls=:dash,lc=:gray10,lw=2)
    vline!([c_crit2],ls=:dash,lc=:blue,lw=2)
    ylabel!(pl,L"\mathrm{max}(x_1)")
    plot!(xlabel=L"w")
    display(pl)
end

savefig(pl,"../../paper_figs/continuation_around_hopf_multiple_attractors_N5_g$(g).png")

#=
params = contres_fp_sols[:,1]
xs = contres_fp_sols[:,2]
@show extrema(xs)
stability = Bool.(contres_fp_sols[:,end])
@show extrema(stability)

linestyle = [stb ? :solid : :dash for stb in stability]
linewidth = [stb ? 1 : 1 for stb in stability]

#-------------------plot fixed point results----------------------
pl = plot(params,xs; lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=linestyle,lw=linewidth,
    xformatter=:auto,widen=true,
    DEFAULTS...,tickfontsize=tick_fontsize+3,yformatter=:auto)

#display(pl)
@show length(contres_fp_bifp)
#plot some (not all) bif points
for bp in eachrow(contres_fp_bifp)
    type,crit_param,idx = bp
    if type == "bp"  
        scatter!([crit_param], [xs[idx]],ms=ms_bp,mc=:orange,markerstrokewidth=0.0)
    elseif type == "nd"
        scatter!([crit_param], [xs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
    end
end


=#