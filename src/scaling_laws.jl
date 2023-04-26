function period_scaling(ds,params;T,Ttr,Δt,state,diffeq,ϵ)
	times = 0:Δt:T
	periods = zeros(length(params))
	new_params = get_parameters(ds)
	
	for (i,p) in enumerate(params)
		@show i,p
		new_params[2] = p
		ts = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state=state,p=new_params,diffeq=diffeq)
		periods[i] = estimate_period(ts[:,1],:autocorrelation,times;L=length(ts[:,1])÷2,ϵ=ϵ)
	end
	return periods
end


function amplitude_scaling(ds,params;T,Ttr,Δt,state,diffeq)
	amplitudes = zeros(length(params))
	new_params = get_parameters(ds)
	
	for (i,p) in enumerate(params)
		@show i,p
		new_params[2] = p
		ts = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state=state,p=new_params,diffeq=diffeq)
		amplitudes[i] = (maximum(ts[:,1]) - minimum(ts[:,1]))/2
	end
	return amplitudes
end
