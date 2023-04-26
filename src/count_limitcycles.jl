function count_limitcycles(ds,ensemble;T,Δt,threshold,Ttr,diffeq)
	#detect period
	#encode limit cycle (unique order?)
	#only powers of 2
	#check if the power sequence is unique
	#1345927689
	N = Int(1/2*dimension(ds))
	sequences = []
	
	for i in 1:ensemble
		state = rand(dimension(ds))
		integer_ts = map_system(ds,T,state;Δt,threshold,Ttr,diffeq)
		sequence,isperiodic = find_limitcycle(integer_ts,N)
		
		if isperiodic == false
			@warn "Timeseries is not periodic! Aborting"
			break
		end
		
		sequence in sequences || push!(sequences,sequence)
	end
	return length(sequences)
end

function find_limitcycle(integer_timeseries,N)
	#need to look at N periods only
	#look at N integer chunks, check if all the same
	#if not, isperiodic = false
	
	
	integer_timeseries = filter_repetitions(integer_timeseries)
	try 
		integer_timeseries = Int8.(log2.(integer_timeseries))
	catch error
		@warn "Transient or parameter value needs to be changed. Integer timeseries has states other than pure states."
		return NaN,false
	end
	
	l = length(integer_timeseries)
	start_index = 1
	powers = 10 .^ [(N-1):-1:0;]
	isperiodic = true
	seq = 0
	
	while l - start_index > N
		Nslice = integer_timeseries[start_index:start_index+N-1]
		seq = Nslice
		
		if start_index != 1 && Nslice != seq
			isperiodic = false
			return NaN,isperiodic
		end
		start_index += N
	end
	
	return seq,isperiodic
end

function findall_limcycles(ds;max_ensemble = 10000,T,Δt,threshold,Ttr,diffeq)
	limcycle_histogram = Dict()
	sequences = Dict()
	nr_unique_limcycles = []
	nr_limcycles = 0
	
	N = Int(1/2*dimension(ds))
	
	for iter in 1:max_ensemble
		iter % 1000 == 0 && @show iter
		
		state = rand(dimension(ds))
		integer_ts = map_system(ds,T,state;Δt,threshold,Ttr,diffeq)
		sequence,isperiodic = find_limitcycle(integer_ts,N)
		
		if isperiodic == false
			@warn "Timeseries is not periodic! Aborting"
			break
		end
		
		
		cycle_index = search_cycle(sequence,sequences)
		
		if cycle_index == nothing
			index_of_new = length(sequences)+1
			push!(sequences,index_of_new => all_circshifts(sequence))
			push!(limcycle_histogram,index_of_new => 1)
			nr_limcycles += 1
			push!(nr_unique_limcycles, nr_limcycles)
		else
			limcycle_histogram[cycle_index] += 1
			push!(nr_unique_limcycles, nr_limcycles)
		end
		
	end
		
	return nr_unique_limcycles,limcycle_histogram

end


function findall_limcycles_many_systems(system_sizes,g;max_ensemble = 10000,T,Δt,threshold,Ttr,diffeq)

	pl1 = plot()
	pl2 = plot()
		
	for N in system_sizes
		@show N
		wc1 = (N-1)/g
		wc2 = (N-1)/N
		w = (wc1+wc2)/2
		
		p = [g,w,N]
		u0 = rand(2*N)
		ds = create_ds(u0,p)
		
		nr_unique_limcycles,limcycle_histogram = findall_limcycles(ds;max_ensemble,T=T,Δt=Δt,threshold=threshold,Ttr=Ttr,diffeq=diffeq)
		
		plot_limc!(pl1,nr_unique_limcycles,factorial(N-1))
		limc_histogram!(pl2,limcycle_histogram)
	end
		
	savefig(pl1,"nr_limcycles_g_$g"*"_max_ensemble_$max_ensemble"*"_Δt_$Δt"*"_th_$threshold"*"_Ttr_$Ttr"*".png")
	savefig(pl2,"limcycles_histogram_g_$g"*"_max_ensemble_$max_ensemble"*"_Δt_$Δt"*"_th_$threshold"*"_Ttr_$Ttr"*".png")

end

all_circshifts(s) = [circshift(s,i) for i in 1:length(s)] 

function search_cycle(c,seqs)
       for s in keys(seqs)
           c in seqs[s] && (return s)
       end
end




