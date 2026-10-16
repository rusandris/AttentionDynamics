using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
using CSV,DataFrames
cd(@__DIR__)

include("../../scripts/new/load_dataset_param_examples.jl")

src_dir_root = "../../data/lyapcurves_datasets_allparams_T500_res100/"
src_files = readdir(src_dir_root) #filter(f->!(isfile(f)),src_dir_root .* readdir(src_dir_root))

#read all dataset parameters
param_datasets_table_file = readdir("../../data/lyapcurve_datasets_parameters")[1]
parameters_dataset = CSV.read("../../data/lyapcurve_datasets_parameters/" * param_datasets_table_file,DataFrame)
αs = []

#dataset groups 
dataset_groups = [["Books","Twitter","Publications"],["Google","Reddit","Movies"]]

markershapes = [:star,:hexagon, :cross, :xcross] 
colors = [:red,:lightseagreen,:orange,:cyan,:brown]

for dataset_group in dataset_groups

    for dataset_category in dataset_group 
        dataset_files = findall(f -> occursin(dataset_category, f), src_files)
        dataset_p = parameters_dataset[occursin.(dataset_category, parameters_dataset[:, :Dataset]), :]
        α = dataset_p[1,:alpha]
        N = dataset_p[1, :N]
        pl = plot()

        for (i, file) in enumerate(dataset_files)
            @show file 
            dataset_name = dataset_p[i]

            lyap_data = readdlm(src_dir_root * file)
            
            #order and save
            ws = lyap_data[:,3]
            order_vals = sortperm(ws)
            ws .= ws[order_vals] #reorder 
            λs = lyap_data[:,end][order_vals] ./ α #scale
            @show length(λs)

            #plot lyap curves 
            plot!(pl,ws,λs,label="")

            #plot param value 
            #λ_pidx = findfirst(w->w==ws_datasets[i],ws) 
            vline!([ws_datasets[i]],lw=2,ls=:dash,lc=colors[i],label=labels[i])

        end
    end
end

plot!(pl_lyapcurves,dpi=400)
lens!(pl_lyapcurves,[0.0,0.07], [0.0, 35.0], inset = (1, bbox(0.5, 0.35, 0.45, 0.45)))
savefig(pl_lyapcurves,"../../paper_figs/lyapcurves_dataset_param_examples_maxima_T_Ttr500.png")
