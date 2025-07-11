function simple_transient_time(ds,state=current_state(ds); Tmax=1000, threshold, period_tol=1e-3, max_tol=1e-3)

	reinit!(ds, state) #reinitialize integrator
	println("id=$(Threads.threadid()): i got " * "$(state[1])" * " and was set" * "$(current_state(ds)[1])" * "my initial time is " * "$(current_time(ds))" )
	previous_period = 0	      #init previous time between peaks of limcycle
	previous_max_timepoint = 0     #init timepoint of previous peak
	prev_max_of_limit_cycle = 0    #init previous peak
	max_timepoint = 0	      #init timepoint of current peak
	max_of_limit_cycle = threshold #init of current peak
	
	transient = Inf	#init of transient 
	found_transient = false 

	
	while current_time(ds) < Tmax	#stepping to Tmax if transient is not found
		step!(ds)	#step the integrator
        state = current_state(ds)
		if state[1] > threshold	#check if value of one variable is not 0ish	
		
			if state[1] > max_of_limit_cycle      #find max(peak)
				max_of_limit_cycle = state[1]
				max_timepoint = current_time(ds)
	
			elseif  max_timepoint != previous_max_timepoint	#check if the integrator left a peak and the timepoints of peaks haven't been updated 
				period = max_timepoint - previous_max_timepoint  
				period_diff = abs(previous_period - period)
				peak_diff = abs(max_of_limit_cycle - prev_max_of_limit_cycle)
				#@show period_diff
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
	println("id=$(Threads.threadid()) finished : final state is " * "$(current_state(ds)[1])")
	return transient
end

#threaded version
function transient_lifetime_statistics(ds,init_states;
		 Tmax=10000, threshold=1.0, period_tol=1e-3, max_tol=1e-3)
	
	ensemble_size = length(init_states)
	#make copies of ds 
	#equal to the number of threads
	ds_copies = [deepcopy(ds) for i in 1:Threads.nthreads()]

	#set different initial conditions to each
	#for	i in 1:length(ds_copies) 
	#	reinit!(ds_copies,init_states[i])
	#end	

	lifetimes = zeros(ensemble_size)
	#calc lifetimes using multithreaded loop
	Threads.@threads for i in 1:ensemble_size
		local id = Threads.threadid()
		println("id : $(id)" ) #* "state: $(current_state(ds_copies[id]))")

		lifetimes[i] = simple_transient_time(ds_copies[id],init_states[i];Tmax=Tmax, threshold=threshold,
        	period_tol=period_tol,max_tol=max_tol)

		@show current_time(ds_copies[id])
		@show current_parameters(ds_copies[id])
		
	end
	return lifetimes
end