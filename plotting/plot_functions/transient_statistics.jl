function plot_histogram(values;bins=100,filename="transient_time_ditribution.png")
	h = histogram(values,
	bins=bins,
	normalize=:probability,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	xlabel=L"$\tau$",
	ylabel=L"$P(\tau)$",
	dpi=300)
	
	savefig(h,filename)

end


function construct_timeseries(trs;n,Δt)
	trs = filter(x -> x != Inf,vec(trs))    #filter out Infs
	lifetime_ts = [n]
	time = [minimum(trs):Δt:maximum(trs)+Δt;]  #construct equidistant timepoints with Δt between
	
	for i in 1:length(time)-1
		escaped = findall(x -> (x >= time[i] && x < time[i+1]),trs) #binning
		nr_escaped = length(escaped)
		
		if nr_escaped != 0 	#check if some transient time values fall into a bin
			push!(lifetime_ts,lifetime_ts[end] - nr_escaped) #substract the escaped
			deleteat!(trs,escaped) #delete the escaped ones from trs
		else
			push!(lifetime_ts,lifetime_ts[end])
		end
	end

	lifetime_ts,time
end

function plot_ensemble_timeseries(time,lifetime_ts;filename="trs_timeseries.png")
	filter!(x -> x != 0,lifetime_ts)
	time = time[1:length(lifetime_ts)]
	
	pl = plot(time,lifetime_ts,
	lc=:gray10,
	lw=2,
	yaxis=:log,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	ylabel=L"$N(t)$",
	xlabel=L"$t$",
	label="",
	ylims=[1,lifetime_ts[1]],
	yticks = [1,10,100,1000],
	dpi=300)	

	a,b = my_linreg(time,log.(lifetime_ts))
	avg_lifetime = round(-1/a,digits=2)
	
	
	plot!(time,exp.(a*time .+ b),
	ls=:dash,
	lc=:red,
	lw=1,
	la=0.5,
	label="$avg_lifetime")
	
	savefig(pl,filename)

end

function plot_ensemble_timeseries(filename;Δt=1)
	trs = readdlm(filename)
	n = length(trs)
	trs_timeseries,time = construct_timeseries(trs;n,Δt)
	plot_ensemble_timeseries(time,trs_timeseries;filename=filename * ".png")
end


function plot_transient_scaling(parameter_values,n,Δt;axistype=:lin,filename="avg_transient_scaling.png")
	avg_lifetimes = []
	
	for f in readdir()
		trs = readdlm(f)
		lifetime_ts,time = construct_timeseries(trs;n=n,Δt=Δt)
		filter!(x -> x != 0,lifetime_ts)
		time = time[1:length(lifetime_ts)]
		
		a,b = my_linreg(time,log.(lifetime_ts))
		avg_lifetime = round(-1/a,digits=2)
		
		push!(avg_lifetimes,avg_lifetime)
		
	end
	
	
	pl = plot(parameter_values,	avg_lifetimes,
	st=:scatter,
	mc=:gray10,
	ms=2,
	xguidefontsize=18,
	tickfontsize=10,
	yguidefontsize=18,
	ylabel=L"$<\tau>$",
	xlabel=L"$w$",
	label="",
	xticks=parameter_values[1:10:end],
	dpi=300)	
	
	a,b = my_linreg(log.(parameter_values),log.(avg_lifetimes))
	@show a,b
	exponent = round(a,digits=2)

	plot!(parameter_values,exp(b)*parameter_values .^ a,
	ls=:dash,
	lc=:red,
	lw=1,
	la=0.5,
	label="$exponent")
	
	if axistype==:log
		plot!(xaxis=:log,yaxis=:log)
	end
	
	savefig(pl,filename)

end

function my_linreg(X,Y)
	x = y = xx = xy = 0
	n=length(X)
	for i in 1:n
		x += X[i]
		xx += X[i]^2
		y += Y[i]
		xy += Y[i]*X[i]
	end
	a=(x*y-n*xy)/(x^2-n*xx)
	b=(x*xy-y*xx)/(x^2-n*xx)

	return (a,b)
end



