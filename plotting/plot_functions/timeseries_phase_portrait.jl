function plot_timeseries!(pl,ts;kwargs...)
	keyword_args = NamedTuple(kwargs)
	T = keyword_args.T
	Ttr = keyword_args.Ttr
	Δt = keyword_args.Δt
	
	time = Ttr:Δt:Ttr+T
	timeseries = plot(time,ts[:,1:2:end];kwargs...)
end

plot_timeseries(ts;kwargs...) = plot_timeseries!(plot(),ts;kwargs...)

function plot_phase_portrait!(pl,xts,yts;kwargs...)
	plot!(pl,xts,yts;kwargs...)
end

plot_phase_portrait(xts,yts;kwargs...) = plot_phase_portrait!(plot(),xts,yts;kwargs...)



