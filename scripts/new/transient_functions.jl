function simple_transient_time(ds,state=current_state(ds); Tmax=1000, threshold, period_tol=1e-3, max_tol=1e-3, return_ts = false)
	reinit!(ds, state) #reinitialize integrator
	previous_period = 0	      #init previous time between peaks of limcycle
	previous_max_timepoint = 0     #init timepoint of previous peak
	prev_max_of_limit_cycle = 0    #init previous peak
	max_timepoint = 0	      #init timepoint of current peak
	max_of_limit_cycle = threshold #init of current peak
	
	transient = Inf	#init of transient 
	found_transient = false 
	
	if return_ts	#return timeseries option
		ts = []
		time = []
	end
	
	while current_time(ds) < T	#stepping to T if transient is not found
		step!(ds)	#step the integrator
		#=
		if return_ts     #save timeseries values if return_ts is true 
			push!(ts,integ.u[1])
			push!(time,integ.t)
		end
		=#
        state = current_state(ds)
		if state[1] > threshold	#check if value of one variable is not 0ish	
		
			if state[1] > max_of_limit_cycle      #find max(peak)
				max_of_limit_cycle = state[1]
				max_timepoint = current_time(ds)
	
			elseif  max_timepoint != previous_max_timepoint	#check if the integrator left a peak and the timepoints of peaks haven't been updated 
				period = max_timepoint - previous_max_timepoint  
				period_diff = abs(previous_period - period)
				peak_diff = abs(max_of_limit_cycle - prev_max_of_limit_cycle)
				@show period_diff
				#@show peak_diff
				#@show previous_period
				#@show period
				#@show prev_max_of_limit_cycle
				#@show max_of_limit_cycle 
				if period_diff < period_tol && peak_diff < max_tol #check for convergence of peak value and period
					transient = current_time(ds)
					found_transient = true
				end
				
				#update previous variables
				previous_period = period
				previous_max_timepoint = max_timepoint
				prev_max_of_limit_cycle = max_of_limit_cycle
			end
		else
			max_of_limit_cycle = threshold   #reset before the integrator arrives at the next peak
		end
		
		found_transient && break      #break if transient is found
	
	end
	
	if return_ts	#return options
		return transient,time,ts
	else 
		return transient
	end
end