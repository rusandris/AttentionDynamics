using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("..//defaults.jl")

cont_dir = "../../data/bifkit_data_g20_N5_trap/"

#create dir for figs
figs_dir = "../paper_figs/poster_figs/"
mkpath(figs_dir)

#saved params
global_params = readdlm(cont_dir*"global_params_Ns_g.txt")
Ns = Int.(global_params[1:end-1])
gval_phase_diag = global_params[end]

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

#read cont results
contres_files = readdir(cont_dir)
branch_plots = []

#find files that match N 
contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]

#fixed point continuation results
contres_fp_sol_file = contres_files_N[findfirst(f -> occursin("contres_fp_sols",f),contres_files_N)]
contres_fp_bifp_file = contres_files_N[findfirst(f -> occursin("contres_fp_bifp",f),contres_files_N)]

#read the matched files
contres_fp_sols = readdlm(cont_dir*contres_fp_sol_file)
contres_fp_bifp = readdlm(cont_dir*contres_fp_bifp_file)

#limit cycle continuation results
contres_sol_files = contres_files_N[findall(f -> occursin("contres_po_sols",f),contres_files_N)]
contres_bifp_files = contres_files_N[findall(f -> occursin("contres_po_bifp",f),contres_files_N)]

#read the matched files
contres_sols = readdlm.(cont_dir .* contres_sol_files)
contres_bifp = readdlm.(cont_dir .* contres_bifp_files)

#critical values
N = 5
g = 20
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
pl = plot(; lc=:gray10,legend=false,guidefontsize=guide_fontsize,
    ls=linestyle,lw=linewidth,
    xformatter=:auto,widen=true,
    DEFAULTS...,tickfontsize=tick_fontsize+3,yformatter=:auto)

#display(pl)
@show length(contres_fp_bifp)
#plot some (not all) bif points
for bp in eachrow(contres_fp_bifp)
    type,crit_param,idx = bp
    if type == "ns"  
        scatter!([crit_param], [xs[idx]],ms=ms_bp,mc=:orange,markerstrokewidth=0.0)
    elseif type == "nd"
        #scatter!([crit_param], [xs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
    end
end

#-------------------plot lim cycle results----------------------
#loop through the two parts of the branch

for i in 1:2
    limc_sols = contres_sols[i]
    limc_bifp = contres_bifp[i]

    params = limc_sols[:,1]
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

end


#display(pl)
#plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

vline!([wc1],ls=:dash,lc=:gray10,lw=2)
vline!([wc2],ls=:dash,lc=:blue,lw=2)
ylabel!(pl,L"\mathrm{max}(x_1)")
display(pl)

#=
#annotations
if i == 1 
    annotate!(0.3, 4.0 , text("limit cycles", :gray10, :center, annotation_size))
elseif i == 2
    annotate!(0.35, 8.0 , text("limit \n cycles", :gray10, :center, annotation_size))
    annotate!(0.25, 3.0 , text("tori", :gray10, :center, annotation_size))
    annotate!(0.25, 12.0 , text("NS", :gray10, :center, annotation_size))

    #annotate!(0.25, 12.0 , text("shorter chaotic \n transients", :gray10, :center, annotation_size))
elseif i == 3
    annotate!(0.32, 29.0 , text("NS", :gray10, :center, annotation_size))
    annotate!(0.28, 26.0 , text("NS", :gray10, :center, annotation_size))
    annotate!(0.22, 19.0 , text("NS", :gray10, :center, annotation_size))

    annotate!(0.38, 20.0 , text("limit \n cycles", :gray10, :center, annotation_size))
    annotate!(0.28, 10.0 , text("tori", :gray10, :center, annotation_size))
    #annotate!(0.25, 24.0 , text("longer chaotic \n transients", :gray10, :center, annotation_size))
end
=#