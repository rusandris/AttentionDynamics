export precompute_trajectories,save_trajectory

function precompute_trajectories(ds, nr_ics, ps::Vector{Float64}, p_idx, Δt=1e-1, T_traj=5000, Ttr=1000)

    for p in ps

        if p_idx isa Vector
            for i in p_idx
                set_parameter!(ds, i, p)
            end
        else
            set_parameter!(ds, p_idx, p)
        end

        #solve systems to get trajectory and save
        precompute_trajectories(ds, nr_ics; T_traj=T_traj, Ttr=Ttr)
    end
end


function precompute_trajectories(ds, nr_ics, ps::Vector{Vector{Float64}}, p_idx::Vector{Int64}; Δt=1e-1,
    T_traj=5000, Ttr=1000, rand_ic=dim -> rand(dim) * 0.1, savedir="./", p_names = ["g","w"],
    filename="traj_test.dat", phasediag_ids=[], ic_counter=0)
    add_ids = isempty(phasediag_ids) ? false : true 

    #loop through parameter pairs
    for (idx,p) in enumerate(ps)
        #@show p
        for (i,idx) in enumerate(p_idx)
            set_parameter!(ds, idx, p[i])
        end

        #put strings together
        par_names = ""
        for i in 1:length(p_names)
            par_names *= p_names[i] * "$(p[i])_"
        end 

        filename_edited = "trajs_" 
        if add_ids 
            id = phasediag_ids[idx]
            filename_edited *= id   
        end
        filename_edited *= "_" * par_names * filename

        #solve systems to get trajectory and save
        precompute_trajectories(ds, nr_ics; T_traj=T_traj, Ttr=Ttr, Δt=Δt, savedir=savedir,
            filename = filename_edited, rand_ic=rand_ic, ic_counter=ic_counter)
    end
end

function precompute_trajectories(ds, nr_ics; savedir="./", filename="traj_test.dat", T_traj=T_traj, ic_counter=0 , Ttr=Ttr, Δt=1e-1, rand_ic=dim -> rand(dim) * 0.1)

    dim = dimension(ds)
    u0s = [rand_ic(dim) for i in 1:nr_ics]
    #@show current_parameters(ds)
    for (i, u0) in enumerate(u0s)
        @info "Initial condition $(i)..."
        traj, times = trajectory(ds, T_traj, u0; Ttr=Ttr, Δt=Δt)
        fname, ext = splitext(filename)
        save_trajectory(traj, fname * "_ic$(i+ic_counter)" * "$(ext)", savedir)
    end

end

function save_trajectory(traj, filename, savedir)
    dim = length(traj[1])
    traj_mat = Matrix(traj[:, 1:2:dim])
    mkpath(savedir)
    write(savedir * filename, vcat(traj_mat...))
end























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

