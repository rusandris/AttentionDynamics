function plot_bifdiagram!(pl,continuation_datafile::String,bifpoints_datafile::String;kwargs...)
	continuation_data = readdlm(continuation_datafile)
	continuation_params = continuation_data[:,1] 
	continuation_curve = continuation_data[:,2]
	
	bifpoints_data = readdlm(bifpoints_datafile)
	bifpoints = [(type = Symbol(bp[1]),param = bp[2],stability_change = Symbol(bp[3])) for bp in eachrow(bifpoints_data)]
	
	#plot branch
	start_index = 1
	for bp in bifpoints
		@show bp
		
		#--------------------plot branch up to bifpoint--------------	
		index_before_bif = findlast(x -> x < bp.param,continuation_params)
		
		#unstable to stable 
		if bp.stability_change == :UnstableToStable ||  bp.stability_change == :UnstableToUnstable 
			ls = :dash
		#stable to unstable
		elseif bp.stability_change == :StableToUnstable
			ls = :solid
		end
		println("Plotting branch...")
		plot!(pl,continuation_params[start_index:index_before_bif],continuation_curve[start_index:index_before_bif];kwargs...,ls=ls)
		
		#--------------------plot bifpoint-------------------------
		bp_color = BIF_COLORS[bp.type]
		bp_coords = [bp.param,continuation_curve[index_before_bif+1]]
		scatter!(pl,[bp_coords[1]],[bp_coords[2]];kwargs...,mc=bp_color)
		annotate!(bp_coords[1]+0.02,bp_coords[2]-0.3,text(BIF_NAMES[bp.type], 12,:left))
	
		start_index = index_before_bif +  2
	end
		ls = bifpoints[end].stability_change == :UnstableToStable ? :solid : :dash
		plot!(pl,continuation_params[start_index:end],continuation_curve[start_index:end];kwargs...,ls=ls)

	return pl
end

plot_bifdiagram(continuation_datafile::String,bifpoints_datafile::String;kwargs...) = 
plot_bifdiagram!(plot(),continuation_datafile::String,bifpoints_datafile::String;kwargs...)

function plot_bifdiagram!(pl,continuation_datafile::String,stability::Symbol;kwargs...)
	continuation_data = readdlm(continuation_datafile)
	continuation_params = continuation_data[:,1] 
	continuation_curve = continuation_data[:,2]
	
	if stability == :Unstable
		ls = :dash
	elseif stability == :Stable
		ls = :solid
	end
	
	plot!(pl,continuation_params,continuation_curve;kwargs...,ls=ls)

end

plot_bifdiagram(continuation_datafile::String,stability::Symbol;kwargs...) = 
plot_bifdiagram!(plot(),continuation_datafile::String,stability::Symbol;kwargs...)
