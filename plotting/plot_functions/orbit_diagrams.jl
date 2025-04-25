function plot_orbit_diagram(filename::String,pvalues;rev=false,ms,ma,saveto)
	od_data = readdlm(filename)
	od_data = [filter(x -> x != "",od_data[i,:]) for i in 1:size(od_data)[1]]
	
	rev && reverse!(od_data,dims=1)
		
	parameters = [fill(pvalues[i],length(od_data[i])) for i in 1:length(pvalues)]
	
	pl = plot(parameters,od_data,
	st=:scatter,
	mc=:gray10,
	legend=false,
	ms=ms,
	markerstrokewidth=0.001,
	markeralpha=ma,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	xlabel=L"$w$",
	ylabel=L"$y_1$",
	dpi=300);

	savefig(pl,saveto)
end


function plot_orbit_diagram(pvalues,poincare_points;kwargs...)
	pl = plot()
	for (i,pval) in enumerate(pvalues)
		points_i = poincare_points[:,i]
		plot!(pl,fill(pval,length(points_i)),points_i;kwargs...)
	end
	return pl
end


function plot_all_orbit_diagrams(dirname, pvalues;saveto)
	cd(dirname)
	filenames = readdir()
	
	pl = plot()
	for f in filenames
		od_data = readdlm(f)
		od_data = [filter(x -> x != "",od_data[i,:]) for i in 1:size(od_data)[1]]
		parameters = [fill(pvalues[i],length(od_data[i])) for i in 1:length(pvalues)]
		
		plot!(parameters,od_data,
		st=:scatter,
		legend=false,
		ms=0.4,
		markerstrokewidth=0.001,
		markeralpha=0.02,
		xguidefontsize=18,
		yguidefontsize=18,
		tickfontsize=10,
		xlabel=L"$w$",
		ylabel=L"$y_1$",
		dpi=400);
	end
	
	savefig(pl,saveto)

end

function plot_2dir_od(filename,pvalues;ms,ma)
	fn,ext = splitext(filename)
	nr_pvalues = length(pvalues)
	od_data = readdlm(filename)
	od_data = [filter(x -> x != "",od_data[i,:]) for i in 1:size(od_data)[1]]

	increasing_od_data = od_data[1:nr_pvalues]
	decreasing_od_data = od_data[end:-1:nr_pvalues+1]
	
	pl = plot()
	
	plot_orbit_diagram!(pl,increasing_od_data,pvalues;ms=ms,ma=ma,mc=:gray10)
	plot_orbit_diagram!(pl,decreasing_od_data,pvalues;ms=ms,ma=ma,mc=:red)
	
	savefig(pl,fn * ".png")

end

function plot_integer_orbit_diagram(filename,pvalues;saveto)
	od = readdlm(filename)
	parameters = [fill(p,size(od)[2]) for p in pvalues]
	
	pl = plot(parameters,collect.(eachrow(od)),
	st=:scatter,
	mc=:gray10,
	legend=false,
	ms=1,
	markerstrokewidth=0.001,
	markeralpha=0.01,
	dpi = 300)
	
	#display(pl)
	savefig(pl,saveto)
end


