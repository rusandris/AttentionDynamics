function period_scaling(ds,params;T,Ttr,Δt,u0=current_state(ds),ϵ)
	periods = zeros(length(params))
	
	for (i,p) in enumerate(params)
		@show i,p
		set_parameter!(ds,2,p)
        traj,t = trajectory(ds,T,u0;Ttr=Ttr)
        x1 = traj[:,1]
		periods[i] = estimate_period(x1,:autocorrelation,t;L=length(x1)÷2,ϵ=ϵ)
	end
	return periods
end

function amplitude_scaling(ds,params;T,Ttr,Δt,u0=current_state(ds))
	amplitudes = zeros(length(params))
	
	for (i,p) in enumerate(params)
		@show i,p
        set_parameter!(ds, 2, p)
        traj, t = trajectory(ds, T, u0; Ttr=Ttr)
        x1 = traj[:, 1]
		amplitudes[i] = (maximum(x1) - minimum(x1))/2
	end
	return amplitudes
end