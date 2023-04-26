function orbit_diagram(ds,plane, idxs, p_index, pvalues; Ttr, tf ,direction = 1,diffeq,rootkw=(xrtol=1e-6,atol=1e-6),filename)
	dim = dimension(ds)

	od = produce_orbitdiagram(ds, plane, idxs, p_index, pvalues; 
		printparams = true,
		direction=direction, 
		warning=true, 
		Ttr=Ttr,
		tfinal = tf, 
		rootkw = rootkw, 
		diffeq=diffeq)
		
	return od
end

function save_od_data(filename, od_data)
	nr_vars = length(od_data[1][1])
	fn,ext = splitext(filename)
	
	for i in 1:nr_vars
		od_data_one_var = []
		for p in 1:length(od_data)
			push!(od_data_one_var,hcat(od_data[p]...)[i,:])
		end
		writedlm(fn * "vars$i"*".txt",od_data_one_var)
	end
end
