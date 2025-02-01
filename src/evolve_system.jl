function evolve_system(ds,T;Ttr,Δt,state,p)
	set_parameters!(ds,p)
	traj,ts = trajectory(ds,T,state;Ttr = Ttr,Δt = Δt)
	
	return traj,ts
end

function evolve_system(ds,t_start,t_stop;Δt,state,p)
	set_parameters!(ds,p)
	traj,time = trajectory(ds,t_stop-t_start,state;Ttr = t_start,Δt = Δt)
	
	return time,traj
end

function ts_to_file(ds,T;Ttr,Δt,state,p)
	w = p[2]
	ts = evolve_system(ds,T;Ttr,Δt,state,p)
	filename = "ts_w_$w"*"_T_$T"*"_Ttr_$Ttr"*"_Δt_$Δt"
	writedlm(filename,ts)
	return filename
end

