function plot_lyapunovs(system_sizes,param_interval)
	pl = plot()

	for (i,f) in enumerate(readdir())
    	lyaps = readdlm(f)
    	size = system_sizes[i]
    	
    	plot!(pl,param_interval,lyaps,
    	label = "N = $size",
    	lw = 2,
    	xlabel = L"$w$",
    	ylabel = L"$\lambda$",
    	xguidefontsize=18,
    	yguidefontsize=18,
    	tickfontsize=10,
    	dpi=300)
    end
    pl
end


function plot_lyapunovs(system_sizes)
	pl = plot()
	colors = distinguishable_colors(length(system_sizes))
	
	for (i,f) in enumerate(readdir())
		lyaps_data = readdlm(f)
		param_interval = lyaps_data[:,1]
		lyaps = lyaps_data[:,2]
		N = system_sizes[i]
		
		diff_from_critical = param_interval .- param_interval[1] 
		
		plot!(pl,diff_from_critical,lyaps,
    	label = "N = $N",
    	lw = 2,
    	lc = colors[i],
    	xlabel = L"$\Delta w$",
    	ylabel = L"$\lambda$",
    	xguidefontsize=18,
    	yguidefontsize=18,
    	tickfontsize=10,
    	dpi=300)
    	
    	wc2_diff = param_interval[end] - param_interval[1]
    	
    	plot!([wc2_diff,wc2_diff],[0,maximum(lyaps)],lw=1,ls=:dash,lc = colors[i],label="")
    	
    	
	
	end
	pl
end

critical_curve(w,N) = (N-1)/w/N

function plot_phase_diagram(Λ,wvalues,gvalues,N;kwargs...)
	pl = plot()
	gvalues = gvalues ./ N

	#plot lyapunov matrix as heatmap
	pl = heatmap(wvalues,gvalues,Λ;kwargs...)
	
	#plot bifurcation curves
	plot!(pl,wvalues,critical_curve.(wvalues,N);kwargs...,ls=:solid)
	plot!(pl,[(N-1)/N,(N-1)/N],[gvalues[1],gvalues[end]];kwargs...,ls=:dot)	
	plot!(ylims=(gvalues[1],gvalues[end]),xlims=(wvalues[1],wvalues[end]))
	return pl
end
