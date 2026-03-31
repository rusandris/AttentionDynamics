using Plots
using LaTeXStrings
using DelimitedFiles
cd(@__DIR__)
include("../plot_functions/orbit_diagrams.jl")
include("defaults.jl")
#-----------------------plot od data-----------------------

#data dirs
#od_dir = "../../orbit_diagram_data_gN5_210/"
#cont_dir = "../../bifkit_data_gN5_210/"
#cont_dir = "../../bifkit_data_gN5/"

od_dir = "../../data/orbit_diagram_data_gN5/"
cont_dir = "../../data/bifkit_data_g5/"

#create dir for figs
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

#saved params
global_params = readdlm(cont_dir*"global_params_Ns_g.txt")
Ns = Int.(global_params[1:end-1])
gval_phase_diag = global_params[end]

#xlims=(0.5,0.65)
xlims_pl=(0.18,0.4)
plot_srate = 1
od_lyap_files = readdir(od_dir)
plots_od = []
tickfontsize_od_bif = 10
yticks_od = [[0:0.5:1;],[0:0.5:1;],[0:5:15;]]
ylims_od = [(0,1.1),(0,1.2),(0,15)]
reduced_margin = -9Plots.mm
reduced_margin_side = -2Plots.mm
ms_bp = 5 
annotation_size = 12

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

    pl = plot_orbit_diagram(ps,od;OD_DEFAULTS...,ma=0.6,ms=0.8,
        ylabel="",xlims=xlims_pl,widen=true,
        xformatter=:none,xlabel="",ylims=ylims_od[i],
        yticks=yticks_od[i],
        tickfontsize=tick_fontsize+3,yformatter=:auto,
        title=L"N=%$N")
        #*"\n"*L"\epsilon(N-1)=%$(round(gval_phase_diag;digits=0))")
    i == 1 && ylabel!(pl,L"x_1")

    vline!([wc1],ls=:dash,lc=:gray10,lw=2)
    vline!([wc2],ls=:dash,lc=:blue,lw=2)

    push!(plots_od,pl)
end
pl_od=plot(plots_od...,layout=(1,length(Ns)));

#----------------plot each branch of periodic orbits---------------------

contres_files = readdir(cont_dir)
yticks_bif = [[0.0:2.0:8.0;],[0.0:6.0:12.0;],[0:5:25;]]
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
    @show g,gval_phase_diag,wc1,wc2

    #--------------------plot bifdiagram------------------

    params = contres_fp_sols[:,1]
    xs = contres_fp_sols[:,2]
    @show extrema(xs)
    stability = Bool.(contres_fp_sols[:,end])
    @show extrema(stability)

    linestyle = [stb ? :solid : :dash for stb in stability]
    linewidth = [stb ? 1 : 1 for stb in stability]

    #-------------------plot fixed point results----------------------
    pl = plot(params, xs; lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        ls=linestyle,lw=linewidth,
        xformatter=:none,xlims=xlims_pl,yticks=yticks_bif[i],widen=true,
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

    params = contres_sols[:,1]
    maxs = contres_sols[:,3]
    mins = contres_sols[:,2]
    stability = Bool.(contres_sols[:,end])

    linestyle = [stb ? :solid : :dot for stb in stability]
    linewidth = [stb ? 4 : 1 for stb in stability]
    #sol_sampling = [stb ? 1 : 20 for stb in stability]

    plot!(pl,params[1:plot_srate:end], maxs[1:plot_srate:end],
        lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        ls=linestyle[1:plot_srate:end],lw=linewidth[1:plot_srate:end])

    #=

    #plot stable and unstable parts separately
    maxs_stable = maxs[stability]  #select stable 
    maxs_unstable = maxs[.!(stability)] #negation

    plot!(pl,params[stability][1:plot_srate:end], maxs_stable[1:plot_srate:end], lc=:gray10,legend=false,guidefontsize=guide_fontsize,
        ls=:solid,lw=1,xformatter=:none,xlims=xlims_pl)

    plot!(pl,params[.!(stability)][1:plot_srate:end], maxs_unstable[1:plot_srate:end], mc=:gray10,legend=false,guidefontsize=guide_fontsize,
        st=:scatter,ms=1,xformatter=:none,xlims=xlims_pl)

    =#

    #plot some (not all) bif points

    for bp in eachrow(contres_bifp)
        type,crit_param,idx = bp
        if type == "ns"  
            scatter!([crit_param], [maxs[idx]],ms=ms_bp,mc=:orange,markerstrokewidth=0.0)
        elseif type == "nd"
            #scatter!([crit_param], [maxs[idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
        elseif type == "bp"
            #scatter!([crit_param], [maxs[idx]],ms=3,mc=:green,markerstrokewidth=0.0)
        end
    end
    #display(pl)
    #plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

    vline!([wc1],ls=:dash,lc=:gray10,lw=2)
    vline!([wc2],ls=:dash,lc=:blue,lw=2)
    i == 1 && ylabel!(pl,L"\mathrm{max}(x_1)")

    #annotations
    if i == 1 
        #annotate!(0.3, 4.0 , text("limit cycles", :gray10, :center, annotation_size))
    elseif i == 2
        #annotate!(0.35, 8.0 , text("limit \n cycles", :gray10, :center, annotation_size))
        #annotate!(0.25, 3.0 , text("tori", :gray10, :center, annotation_size))
        annotate!(0.25, 12.0 , text("NS", :gray10, :center, annotation_size))

        #annotate!(0.25, 12.0 , text("shorter chaotic \n transients", :gray10, :center, annotation_size))
    elseif i == 3
        annotate!(0.32, 29.0 , text("NS", :gray10, :center, annotation_size))
        annotate!(0.28, 26.0 , text("NS", :gray10, :center, annotation_size))
        annotate!(0.22, 19.0 , text("NS", :gray10, :center, annotation_size))

        #annotate!(0.38, 20.0 , text("limit \n cycles", :gray10, :center, annotation_size))
        #annotate!(0.28, 10.0 , text("tori", :gray10, :center, annotation_size))
        #annotate!(0.25, 24.0 , text("longer chaotic \n transients", :gray10, :center, annotation_size))
    end
    push!(branch_plots,pl)
end

#------------------------plot Floquet exponents/multipliers--------------------------
#plot all non-unity floquet multis (real and imag parts)    
#plot trace and det? 
#plot real part?
#plot complex plane with unit circle?
plots_maxfloq = []
ylims_floq = (-1.0,1.0)

for (i,N) in enumerate(Ns)
    #find files that match N 
    contres_files_N = contres_files[findall(f -> occursin("N$N",f),contres_files)]
    contres_eig_file = contres_files_N[findfirst(f -> occursin("contres_po_eig",f),contres_files_N)]
    contres_sol_file = contres_files_N[findfirst(f -> occursin("contres_po_sols",f),contres_files_N)]

    #find files that match N 
    od_lyap_files_N = od_lyap_files[findall(f -> occursin("N$N",f),od_lyap_files)]
    lyaps_file = od_lyap_files_N[findfirst(f -> occursin("lyaps",f),od_lyap_files_N)]

    @show contres_eig_file
    #read the matched files
    contres_eig = readdlm(cont_dir*contres_eig_file)

    #filter eig results
    #contains real and imag parts separately
    #first N cols are real parts
    #last N cols are imag parts

    realparts = contres_eig[:,1:2N]
    imagparts = contres_eig[:,(2N+1):end]

    #how to filter out zero exponents
    #only select floquet multipliers that cross zero 
    extremas = vec(extrema(realparts;dims=1))
    signs = [extremas[i][1]*extremas[i][2] for i in 1:2N]
    realparts_sign_change = signs .< 0.0
    @show realparts_sign_change


    if N > 3
        realparts = realparts[:,realparts_sign_change]
    end
    #realparts_max = maximum(realparts;dims=2)

    @show cont_dir*contres_sol_file
    contres_sols = readdlm(cont_dir*contres_sol_file)
    lyaps = readdlm(od_dir*lyaps_file)[:,2]
    lyaps_params = readdlm(od_dir*lyaps_file)[:,1]

    params = contres_sols[:,1]
    #@show params

    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #floq maxreal
    #@show length(params),length(realparts_max)
    #@show realparts[:,1]
    labels = []
    pl_maxfloq = plot(params,realparts[:,1];xformatter=:auto,
        guidefontsize=guide_fontsize,xlabel=L"w",widen=true,
        ylims=ylims_floq,xlims=xlims_pl,ms=1, lw=2,
        legend_columns=2,legend=:bottomright,
        lc=:gray50,markerstrokewidth=0.0,label=L"\mu_i",DEFAULTS...,
        yformatter=:auto,tickfontsize=tick_fontsize+3,xticks=[0.3,0.4],
        legendfontsize = 10)

    plot!(params,realparts[:,2:end],label="",lc=:gray50)
    #floq maximag
    #plot!(pl_maxfloq,params,imagpart,mc=:gray50,xlabel=L"w",guidefontsize=guide_fontsize,
    #    xlims=xlims_pl,markershape=:cross,ms=1,lc=:gray50,markerstrokewidth=0.0,legend=false)

    #lyapunov LLE
    plot!(lyaps_params,lyaps,lw=3,mc=:gray10,lc=:gray10,markerstrokewidth=0.0,label=L"\Lambda")

    vline!([wc1],ls=:dash,lc=:gray10,lw=2,label="")
    #vline!([wc2],ls=:dash,lc=:blue,lw=2,label="")

    hline!([0.0],ls=:dash,lc=:gray10,lw=0.5,label="")
    i == 1 && ylabel!(pl_maxfloq,L"\mathrm{Re}(\mu),\Lambda")
    yticks!([-1,0,1])

    #annotations
    if i == 1 
        #annotate!(0.3, 4.0 , text("limit cycles", :gray10, :center, annotation_size))
        annotate!(wc1, -1.3, text(L"w_{\mathrm{c}}", :gray10, :center, 16))
    elseif i == 2
        #annotate!(0.35, 0.3 , text("shorter chaotic \n transients", :gray10, :center, annotation_size-5))
        annotate!(wc1, -1.3, text(L"w_{\mathrm{c}}", :gray10, :center, 16))
    elseif i == 3
        #annotate!(0.37, 0.3 , text("longer chaotic \n transients", :gray10, :center, annotation_size-5))
        annotate!(wc1, -1.3, text(L"w_{\mathrm{c}}", :gray10, :center, 16))
    end

    push!(plots_maxfloq,pl_maxfloq)



end

#l = @layout [a{0.5h}; b{0.25h}; c{0.25h} ]



#=

#plot transient lifetimes 

tr_data_dir = "../../data/transient_lifetimes/"
transient_lifetime_files = readdir(tr_data_dir)
plots_trs = []
ylims_trs = (0,1800)

for (i,N) in enumerate(Ns[end:end])
    trs_file_N = transient_lifetime_files[findall(f -> occursin("N$N",f),transient_lifetime_files)][1]
    @show trs_file_N
    transient_data = readdlm(tr_data_dir*trs_file_N) 
    ws = readdlm(tr_data_dir*"w_values.txt")
    transient_avgs = transient_data[:,1]

    pl = plot(log10.(reverse(ws)),reverse(log10.(transient_avgs));
        xlabel=L"w",ylabel="",xformatter=:auto,
        guidefontsize=guide_fontsize,widen=true,
        ms=1, lw=2,
        legend_columns=2,legend=false,
        lc=:gray10,markerstrokewidth=0.0,DEFAULTS...,
        yformatter=:auto,tickfontsize=tick_fontsize+3)

    #i == 1 && ylabel!(pl,L"\langle \tau \rangle")


    #vline!([0.2],ls=:dash,lc=:gray10,lw=2,label="")
    push!(plots_trs,pl)

end
#pl_trs = plot(plots_trs...,layout=(1,length(Ns)))#bottom_margin=reduced_bottom_margin)
=#

pl = plot(plots_od...,branch_plots...,plots_maxfloq...,layout=(3,3),size=(700,400),
    dpi=300)

#add plot label annotation 
annotate!(pl[1], 0.11, 1.3, text(L"(b)", :gray10, :center, 18))


#pl = plot(pl_od,pl_bif,pl_floq,layout=(3,1),size=(1000,500),dpi=300,left_margin=7Plots.mm,bottom_margin=7Plots.mm)
#savefig(pl,"../../continuation_test_maxreal.png")

#savefig(pl,figs_dir*"continuation_results_bif_floq_lyap_trans_od_Ttr20000_noreinit_g5_N310_test.png")

