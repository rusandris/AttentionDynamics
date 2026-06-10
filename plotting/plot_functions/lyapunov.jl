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

#article params in rescaled param space
w_rescaled_param(r,N;α=0.005,K=1.0,c=2.4) = 1/(1+r/(α*K*c*(N-1)))
g_rescaled_param(r;α=0.005) = r/α

critical_curve(w,N) = (N-1)/w
critical_curve_adim(w) = 1/w
w_g_dependence_curve(w,N;K=1.0,c=2.4) = K*c*(N-1)*(1/w -1)
w_g_dependence_curve_adim(w;K=1.0,c=2.4) = K*c*(1/w -1)

function plot_phase_diagram(Λ, wvalues, gvalues, N; plot_critical_curve=true, 
    critical_curve_colors=[:gray10, :gray80], critical_curve_lines=[:dash,:dash],
    critical_curve_lw=[3,3],
    kwargs...)
	pl = plot()
	#plot lyapunov matrix as heatmap
	pl = heatmap(wvalues,gvalues,Λ;kwargs...)
	
	#plot bifurcation curves
	if plot_critical_curve 
        plot!(pl, wvalues, critical_curve_adim.(wvalues); kwargs..., 
            ls=critical_curve_lines[1],
            lw = critical_curve_lw[1],
            lc=critical_curve_colors[1], label="") #label=L"g_c(w)")
        plot!(pl, [(N - 1) / N-0.001, (N - 1) / N-0.001], [gvalues[1], gvalues[end]]; kwargs..., 
            lc=critical_curve_colors[2], 
            ls=critical_curve_lines[2],
            lw=critical_curve_lw[2],
            label="") #label=L"w_{d}=(N-1)/N")	
    end
	#plot!(ylims=(gvalues[1]/(N-1),gvalues[end]/(N-1)),xlims=(wvalues[1],wvalues[end]))
	return pl
end

function plot_article_params(rs,N;ma=[0.4,0.6,0.8,1.0],ylims=[0,20])
	pl = plot(ws,w_g_dependence_curve_adim.(ws),lc=:green,lw=2,ylabel=L"g/(N-1)",
	xlabel=L"w",label=L"g(w)_{r_c=r_p}",legendfontsize=12,ylims=ylims)
	plot!(pl,ws,critical_curve_adim.(ws),lc=:gray10,lw=2,label=L"g_c(w)")

	for (i,r) in enumerate(rs)
		plot!([w(r,N)],[g(r)/(N-1)],st=:scatter,ms=6,mc=:green,ma=ma[i],label=L"r=%$(r)")
	end

	return pl
end