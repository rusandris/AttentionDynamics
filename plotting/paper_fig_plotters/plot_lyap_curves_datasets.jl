using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
cd(@__DIR__)

include("../../scripts/new/load_dataset_param_examples.jl")

src_dir_root = "../../data/lyapcurves_datasets_T500_twores/"
src_dirs = filter(f->!(isfile(f)),src_dir_root .* readdir(src_dir_root))[[2,3,4,5,1]]
labels = ["Movies 1995-2000","Publications 1990-2010","Reddit 2012","Google Trends 2014","Books 1930-1950"]
αs = [p_const_mov[5],p_const_pub[5],p_const_reddit[5],p_const_trends[5],p_const_books[5]]
ys_datasets = [y_mov,y_pub,y_red,y_tr,y_books]
ws_datasets = [w_mov,w_pub,w_red,w_tr,w_books]
markershapes = [:star,:hexagon, :cross, :xcross] 
colors = [:red,:lightseagreen,:orange,:cyan,:brown]

pl_lyapcurves = plot(guidefontsize=18,
    tickfontsize=15,
    title=L"T/\alpha = T_{tr}/\alpha = 500",
    xlabel=L"w",ylabel=L"\Lambda")

for (i,dir) in enumerate(src_dirs)
    
    #read file in dir
    f = readdir(dir)[1]
    lyap_data = readdlm(dir * "/" * f)
 
    #order and save
    ws = lyap_data[:,3]
    order_vals = sortperm(ws)
    ws .= ws[order_vals] #reorder 
    λs = lyap_data[:,end][order_vals] ./ αs[i] #scale
    @show length(λs)

    #plot lyap curves 
    plot!(pl_lyapcurves,ws,λs,label="",lc=colors[i])

    #plot param value 
    #λ_pidx = findfirst(w->w==ws_datasets[i],ws) 
    vline!([ws_datasets[i]],lw=2,ls=:dash,lc=colors[i],label=labels[i])

end

plot!(pl_lyapcurves,dpi=400)
lens!(pl_lyapcurves,[0.0,0.07], [0.0, 35.0], inset = (1, bbox(0.5, 0.35, 0.45, 0.45)))
savefig(pl_lyapcurves,"../../paper_figs/lyapcurves_dataset_param_examples_maxima_T_Ttr500.png")
