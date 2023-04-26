#turn the 2N dimensional dynamical system into a map
#such that n number of variables exceeding a threshold value (active)
#corresponds to a binary state vector

function map_to_integer(state,threshold)
	binary_state = state .>= threshold
	integer = convert_to_int(binary_state)
end

function map_to_binarystate(state,threshold)
	binary_state = state .>= threshold
end

function convert_to_int(binary_state)
	integer = 0
	for i in 1:length(binary_state)
		binary_state[i] && (integer += 2^(i-1))
	end
	integer
end

function map_timeseries_to_binarystate(ds,T,state;Δt,threshold,Ttr,diffeq,return_integer)
	timeseries = trajectory(ds,T,state;Ttr = Ttr,Δt=Δt,save_idxs=1:2:dimension(ds),diffeq=diffeq)
	binary_ts = map_to_binarystate.(eachrow(timeseries),threshold)
	filter!(state -> !all(x -> x == 0,state),binary_ts)
	
	length(binary_ts) == 0 && @warn "Threshold is too high. All states are zero states."
	if return_integer
		integer_ts = convert_to_int.(binary_ts)
		return binary_ts,integer_ts
	else
		return binary_ts
	end
	
	binary_ts
end

function map_timeseries_to_integer(ds,T,state;Δt,threshold,Ttr,diffeq)
	timeseries = trajectory(ds,T,state;Ttr = Ttr,Δt=Δt,save_idxs=1:2:dimension(ds),diffeq=diffeq)
	integer_ts = map_to_integer.(eachrow(timeseries),threshold)
	filter!(x->x != 0,integer_ts)
	
	length(integer_ts) == 0 && @warn "Threshold is too high. All states are zero states."
	
	integer_ts
end

function generate_randomstates(nr,N)
	state = rand([true,false],N)
	random_states = []
	
	for n in 1:nr
		new_state = rand([true,false],N)
		if new_state == state || all(x->x==0,new_state)
			while new_state == state || all(x->x==0,new_state)
				new_state = rand([true,false],N)
			end
		end
		push!(random_states,convert_to_int(new_state))
		state = new_state
	end
	random_states
end

function flip!(state,flipped)
	positions = sample(1:length(state),flipped,replace=false)
	for pos in positions
		state[pos] = !state[pos]
	end
end

function flip_randomstates(nr,N;flipped,return_binary)
	state = rand([true,false],N)
	int_states = []
	bin_states = []
	
	for n in 1:nr
		flip!(state,flipped)
		
		if return_binary
			push!(bin_states,deepcopy(state))
		end
		push!(int_states,convert_to_int(state))
	end
	if return_binary
		return bin_states,int_states
	else
		return int_states
	end
end

function bifurcation_diagram(ds,p_index,pvalues;T,Ttr,init_state = get_state(ds),Δt,threshold,diffeq)
	its_container = []
	state = init_state
	for p in pvalues
		set_parameter!(ds,p_index,p)
		@show p
		its = map_system(ds,T,state;Δt = Δt,threshold = threshold,Ttr = Ttr,diffeq = diffeq)
		push!(its_container,its)
	end
	its_container
end

function filter_repetitions(ts)
	fts = [ts[1]]
	for i in 1:length(ts)-1
		if ts[i] != ts[i+1] 
			push!(fts,ts[i+1])
		else 
			continue
		end
	end
	fts
end
