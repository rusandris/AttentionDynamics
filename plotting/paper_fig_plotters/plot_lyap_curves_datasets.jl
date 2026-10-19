using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
using CSV,DataFrames
cd(@__DIR__)

include("../../scripts/new/load_dataset_param_examples.jl")

src_dir_root = "../../data/lyapcurves_datasets_allparams_T500_res100/"
src_files = readdir(src_dir_root) #filter(f->!(isfile(f)),src_dir_root .* readdir(src_dir_root))

#read all dataset parameters
param_datasets_table_file = readdir("../../data/lyapcurve_datasets_parameters/")[2]
parameters_dataset = CSV.read("../../data/lyapcurve_datasets_parameters/" * param_datasets_table_file,DataFrame)

#dataset groups 
dataset_groups = [["Books","Twitter","Publications"],["Google","Movies","Reddit"]]
lyapcurves_all = []

#read out data
for (gi,dataset_group) in enumerate(dataset_groups)
    
    lyapcurves_group = []

    for (ci,dataset_category) in enumerate(dataset_group)

        dataset_files = src_files[findall(f -> occursin(dataset_category, f), src_files)]
        dataset_p = parameters_dataset[occursin.(dataset_category, parameters_dataset[:, :Dataset]), :]
        α = dataset_p[1,:alpha]
        N = dataset_p[1, :N]
        lyapcurves_cat = []

        for (i, file) in enumerate(dataset_files)
            @show file 

            #read lyapcurves
            lyap_data = readdlm(src_dir_root * file)
            
            #order and save
            ws = lyap_data[:,3]
            order_vals = sortperm(ws)
            ws .= ws[order_vals] #reorder 
            Λs = lyap_data[:,end][order_vals] ./ α #scale
            @show length(Λs),length(ws) 
            push!(lyapcurves_cat,hcat(ws,Λs))

        end
        push!(lyapcurves_group,lyapcurves_cat)
        @show length(lyapcurves_cat)
    end
    push!(lyapcurves_all,lyapcurves_group)
end

#plot lyapcurves and parameters markers separately 
plots = []
mshapes = [[:pentagon,:dtriangle,:circle], [:hexagon, :diamond,:octagon]]
colorpairs = [[[:yellow,:orange],[:darkgreen,:lime],[:darkblue,:lightblue]],
            [[:yellow,:orange],[:darkgreen,:lime],[:darkblue,:lightblue]]]

xlims_groups = [(0.0,0.5),(0.0,0.04)]
ms=6;lw=2

#curves
for (gi,lyapcurves_group) in enumerate(lyapcurves_all)
    @show gi
    pl_group = plot(xlabel=L"w",ylabel=L"\Lambda")
    colors_group = colorpairs[gi]
    
    for (ci,lyapcurves_cat) in enumerate(lyapcurves_group)

        dataset_category = dataset_groups[gi][ci]
        dataset_p = parameters_dataset[occursin.(dataset_category, parameters_dataset[:, :Dataset]), :]
        colors = palette(colors_group[ci],length(lyapcurves_cat)) 
        @show length(colors)
        linealphas = range(0.1,1.0;length=length(lyapcurves_cat))
        msizes = range(6,15;length=length(lyapcurves_cat))

        #plot params afterwards
        for (i, lyapcurve) in enumerate(lyapcurves_cat)

            #plot lyap curves 
            color = colors[i]
            Λs = lyapcurve[:,2]
            ws = lyapcurve[:,1]
            plot!(pl_group,ws,Λs,label="",
                lc=color,
                xlims=xlims_groups[gi],
                lw=lw)

        end
    end
    push!(plots,pl_group)
end

#parameters markers
for (gi,lyapcurves_group) in enumerate(lyapcurves_all)
    colors_group = colorpairs[gi]
    for (ci,lyapcurves_cat) in enumerate(lyapcurves_group)
        dataset_category = dataset_groups[gi][ci]
        dataset_p = parameters_dataset[occursin.(dataset_category, parameters_dataset[:, :Dataset]), :]
        colors = palette(colors_group[ci],length(lyapcurves_cat)) 
        for (i, lyapcurve) in enumerate(lyapcurves_cat)

            #plot param value 
            dataset_name = dataset_p[i,:][1]
            w_dataset = dataset_p[i,:][6]
            ws = lyapcurve[:,1]
            Λs = lyapcurve[:,2]
            #find closest (w,Λ) pair closest to fitparam
            _,closest_w_idx = findmin(abs.( ws .- w_dataset ))

            @show ws[closest_w_idx], w_dataset
            color = colors[i]
            plot!(plots[gi],[ws[closest_w_idx]],[Λs[closest_w_idx]],
                st=scatter,markerstrokewidth=0.5,ms=ms,
                mc=color,markershape=mshapes[gi][ci],label=dataset_name)

            #vline!([w_dataset],ls=:dot,ms=4,
            #    lc=color,label=dataset_name)
        end
    end
end

pl_lyapcurves = plot(plots...,layout=(1,2),size=(1200,400),legend=:outertopright,
    bottom_margin=10Plots.mm,left_margin=8Plots.mm,dpi=600,
    guidefontsize=18,tickfontsize=10)

#lens!(pl_lyapcurves,[0.0,0.07], [0.0, 35.0], inset = (1, bbox(0.5, 0.35, 0.45, 0.45)))
savefig(pl_lyapcurves,"../../paper_figs/lyapcurves_dataset_params_all_maxima_T_Ttr500.png")
