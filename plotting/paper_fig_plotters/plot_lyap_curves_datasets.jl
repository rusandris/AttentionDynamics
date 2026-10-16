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

mshapes = [[:utriangle,:dtriangle,:circle], [:hexagon, :diamond,:octagon]]
colorpairs = [[[:orange,:orangered4],[:lime,:darkgreen],[:lightblue,:darkblue]],[[:orange,:orangered4],[:lime,:darkgreen],[:lightblue,:darkblue]]]
plots = []
xlims_groups = [(0.0,0.5),(0.0,0.04)]
ms=8;lw=2
for (gi,dataset_group) in enumerate(dataset_groups)
    
    colors_group = colorpairs[gi]
    pl_group = plot(xlabel=L"w",ylabel=L"\Lambda")

    for (ci,dataset_category) in enumerate(dataset_group)

        dataset_files = src_files[findall(f -> occursin(dataset_category, f), src_files)]
        dataset_p = parameters_dataset[occursin.(dataset_category, parameters_dataset[:, :Dataset]), :]
        α = dataset_p[1,:alpha]
        N = dataset_p[1, :N]
        colors = palette(colors_group[ci],length(dataset_files)) 
        linealphas = range(0.1,1.0;length=length(dataset_files))
        msizes = range(6,15;length=length(dataset_files))

        for (i, file) in enumerate(dataset_files)
            @show file 
            dataset_name = dataset_p[i,:][1]
            w_dataset = dataset_p[i,:][6]

            #read lyapcurves
            lyap_data = readdlm(src_dir_root * file)
            
            #order and save
            ws = lyap_data[:,3]
            order_vals = sortperm(ws)
            ws .= ws[order_vals] #reorder 
            Λs = lyap_data[:,end][order_vals] ./ α #scale


            #plot lyap curves 
            color = colors[i]
            plot!(pl_group,ws,Λs,label="",lc=color,
                xlims=xlims_groups[gi],
                lw=lw)

            #plot param value 
            #find closest (w,Λ) pair closest to fitparam
            _,closest_w_idx = findmin(abs.( ws .- w_dataset ))

            @show ws[closest_w_idx], w_dataset
            plot!([ws[closest_w_idx]],[Λs[closest_w_idx]],
                st=scatter,ms=msizes,markerstrokewidth=0.5,
                mc=color,markershape=mshapes[gi][ci],label=dataset_name)

            #vline!([w_dataset],ls=:dot,ms=4,
            #    lc=color,label=dataset_name)

        end
    end
    push!(plots,pl_group)
end

pl_lyapcurves = plot(plots...,layout=(1,2),size=(1200,400),legend=:outertopright,
    bottom_margin=10Plots.mm,left_margin=8Plots.mm,dpi=600,
    guidefontsize=18,tickfontsize=10)

#lens!(pl_lyapcurves,[0.0,0.07], [0.0, 35.0], inset = (1, bbox(0.5, 0.35, 0.45, 0.45)))
savefig(pl_lyapcurves,"../../paper_figs/lyapcurves_dataset_params_all_maxima_T_Ttr500.png")
