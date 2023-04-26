function evolve_system(ds,T;Ttr,Δt,state,p,diffeq)
	set_parameter!(ds,p)
	ts = trajectory(ds,T,state;Ttr = Ttr,Δt = Δt,diffeq=diffeq)
	
	return Matrix(ts)
end

function evolve_system(ds,t_start,t_stop;Δt,state,p,diffeq)
	set_parameter!(ds,p)
	time = t_start:Δt:t_stop
	ts = trajectory(ds,t_stop-t_start,state;Ttr = t_start,Δt = Δt,diffeq=diffeq)
	
	return time,Matrix(ts)
end

function ts_to_file(ds,T;Ttr,Δt,state,p,diffeq)
	w = p[2]
	ts = evolve_system(ds,T;Ttr,Δt,state,p,diffeq)
	filename = "ts_w_$w"*"_T_$T"*"_Ttr_$Ttr"*"_Δt_$Δt"
	writedlm(filename,ts)
	return filename
end

