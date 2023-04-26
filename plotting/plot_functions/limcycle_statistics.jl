function plot_limc(nr_limc,max_nr_limcycles,filename="nr_limc.png")
	
	pl = plot(nr_limc,
	lc=:gray10,
	lw = 2,
	ylims=[0.0,1.3*max_nr_limcycles],
	legend=false,
	xguidefontsize=18,
	tickfontsize=8,
	yguidefontsize=18,
	xlabel="Number of iterations",
	ylabel="Number of limit cycles",
	dpi=300)
	
	plot!([0,length(nr_limc)],[max_nr_limcycles,max_nr_limcycles],
	ls=:dash,
	lc=:red,
	lw=2,
	la=1.)
	
	savefig(pl,filename)
end

function plot_limc!(pl,nr_limc,max_nr_limcycles)
	
	plot!(pl,nr_limc,
	legend=false,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	xlabel="Number of iterations",
	ylabel="Number of limit cycle",
	dpi=300)
	
	plot!(pl,[0,length(nr_limc)],[max_nr_limcycles,max_nr_limcycles],
	ls=:dash,
	la=0.5)
end


function limc_histogram(limc_histogram,filename="limc_histogram.png")
	limc_histogram = sort(limc_histogram,by = x -> x[1])
	
	index_limc = collect(keys(limc_histogram))
	freq_limc = collect(values(limc_histogram))
	
	pl = plot(index_limc,freq_limc,
	linestyle = :dashdot,
	lc=:gray10,
	markershape=:circle,
	markercolor=:gray10,
	legend=false,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	ylims = [0,maximum(freq_limc)+5],
	xlabel="Index",
	ylabel="Frequency of limit cycle",
	dpi=300)

	savefig(pl,filename)
end


function limc_histogram!(pl,limc_histogram)
	limc_hist = collect(limc_histogram)
	sort!(limc_hist,by = x -> x[1])
	
	index_limc = [pair[1] for pair in limc_hist]
	freq_limc = [pair[2] for pair in limc_hist]

	plot!(pl,index_limc,freq_limc,
	linestyle = :dashdot,
	markershape=:circle,
	legend=false,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	xlabel="Index",
	ylabel="Frequency of limit cycle",
	dpi=300)

end

