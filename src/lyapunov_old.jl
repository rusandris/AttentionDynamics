export lyapunov_curve,lyapunov_curves,lyapunov_diagram,lyapunov_diagram!
export calc_lyapunov_curves,calc_lyapunov_curves_reinit,ensemble_average_lyapunov

#calculate lyapunov curves for `ws` values of w 
#using several different Ttr transient times 
#for the same g 
#without reinit 
function calc_lyapunov_curves(ds,ws,Ttrs,T;u0=current_state(u0),save_dir="./")

    lyap_curves = [] #save lyaps here

    for Ttr in Ttrs
        g,_,N = current_parameters(ds)
        N = Int(N)
        #create pds 
        pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9 .* randn(2N)])
        file_name = "lyapunovs_g_$g"*"_N_$N"*"_Ttr_$Ttr"*".txt"

        for (i,w) in enumerate(ws)
            println("Ttr: $Ttr" * " w = $w")
            set_parameter!(pds,2,w)
            λ = lyapunov(pds,T;Ttr=Ttr)
            
            open(save_dir*file_name, "a") do f
                write(f, "$w \t $λ \n")
            end

        end
    end
end

#calculate lyapunov curves for `ws` values of w 
#using several different Ttr transient times 
#for the same g 
#with reinit
#paralellized version 
function calc_lyapunov_curves_reinit(ds,ws,Ttrs,T;
        generate_u0 = dim -> generate_initial_condition(dim ÷ 2),d0=1e-9)

    dim = dimension(ds)
    lyap_curves = [] #save lyaps here

    for Ttr in Ttrs
        #save lyapunov exps for a transient time  
        lyaps = zeros(length(ws))

        #distribute parameter loop over threads
        Threads.@threads for i in 1:length(ws)
            
            w = ws[i]
            id = Threads.threadid()
            println("Id: $id" * " Ttr: $Ttr" * " w = $w")

            #make a pds for every param value
            #to avoid races
            #generate new initial condition 
            u0 = generate_u0(dim)
            pds = ParallelDynamicalSystem(deepcopy(ds),[u0,u0 .+ d0 .* rand(dim)])
            set_parameter!(pds,2,w)
            lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
        end
        push!(lyap_curves,lyaps)
    end

    return lyap_curves
end

#initialize ensemble from the chaotic attractor?
#for each parameter, take the average LLE  
function ensemble_average_lyapunov(ds,ws,T,Ttr,nr_ens;u0=current_state(ds),
        ensemble_init_param,T_init_traj = T,Ttr_init_traj = Ttr,d0=1e-9,Δt = 0.01,
        filename="all_lyaps.txt")

    lyap_avgs = zeros(length(ws)) #save lyap avgs here
    dim = dimension(ds)

    #init ensemble from a trajectory
    set_parameter!(ds,2,ensemble_init_param)
    traj,t = trajectory(ds,T_init_traj;Ttr=Ttr_init_traj)
    lyap_exps_all = zeros(length(nr_ens),length(ws))

    for i in 1:length(ws)
    	w = ws[i]
	println("w = $w")
        #go over initial conditions (paralellized)
        #select nr_ens init conds from traj
        lyaps_ensemble = zeros(nr_ens) #save lyap avgs here

        Threads.@threads for e in 1:nr_ens
            id = Threads.threadid()
            #println("Thread Id: $id" * " ensemble member: $e" * " w = $w")

            #initialize from previously computed traj
            ic_idx = rand(1:length(traj))
            u0 = traj[ic_idx]
            pds = ParallelDynamicalSystem(deepcopy(ds),[u0,u0 .+ d0 .* rand(dim)])
            set_parameter!(pds,2,w)
            lyaps_ensemble[e] = lyapunov(pds,T;Ttr=Ttr)
        end
        

	open(filename, "a") do f
	   writedlm(f, lyaps_ensemble')
	end


        lyap_avgs[i] = mean(lyaps_ensemble)
    end

    return lyap_avgs
end

function lyapunov_curve(pinteg,T,Ttr,param_interval;Δt, d0, ut, lt)
	lyaps = []
	for param in param_interval
		reinit!(pinteg)
		set_parameter!(pinteg,2,param)
		lyap = lyapunov(pinteg, T, Ttr, Δt, d0, ut, lt)
		push!(lyaps,lyap)
	end
	return lyaps
end

function lyapunov_curves(system_sizes,param_interval::Vector{Float64};T,Ttr,wc,Δt, d0, ut, lt,diffeq,save=false)
	adapt = diffeq.adaptive
	u0 = rand(2*system_sizes[1]) #dummy initialization
	p = zeros(3)
	ds = create_ds(u0,p)
	
	
	lyap_curves = []
	
	for size in system_sizes	
		u0 = rand(2*size)
		
		pinteg = parallel_integrator(ds,[u0,u0+d0*rand(2*size)];diffeq=diffeq)
		set_parameter!(pinteg,3,size)
		g = (size-1)/wc
		set_parameter!(pinteg,1,g)
		lyaps = lyapunov_curve(pinteg,T,Ttr,param_interval;Δt, d0, ut, lt)
		
		if save 
			filename = "lyaps_samewc_size_$size"*"_T_$T"*"Ttr_$Ttr"*"_adapt_$adapt"*".txt"
			writedlm(filename,lyaps)
			println("creating $filename")
		else
			push!(lyap_curves,lyaps)
		end
	end
	if save == false 
		return lyap_curves
	end
end


function lyapunov_curves(system_sizes,nr_points::Int;T,Ttr,g,Δt, d0, ut, lt,diffeq,save=false)
	adapt = diffeq.adaptive
	u0 = rand(2*system_sizes[1]) #dummy initialization
	p = zeros(3)
	ds = create_ds(u0,p)
	
	
	lyap_curves = []
	
	for N in system_sizes	
		u0 = rand(2*N)
		
		pinteg = parallel_integrator(ds,[u0,u0+d0*rand(2*N)];diffeq=diffeq)
		set_parameter!(pinteg,3,N)
		set_parameter!(pinteg,1,g)
		wc1 = (N-1)/g
		wc2 = (N-1)/N
		param_interval = range(wc1,wc2,nr_points)
		
		lyaps = lyapunov_curve(pinteg,T,Ttr,param_interval;Δt, d0, ut, lt)
		
		if save 
			filename = "lyaps_size_$N"*"_T_$T"*"Ttr_$Ttr"*"_adapt_$adapt"*"_g_$g.txt"
			println("creating $filename")
			writedlm(filename,hcat(param_interval,lyaps))
		else
			push!(lyap_curves,lyaps)
		end
	end
	if save == false 
		return lyap_curves
	end
end

function create_pinteg(u0,p,d0;diffeq)
	ds = create_ds(u0,p)
	N = Int(p[3])
	pinteg = parallel_integrator(ds,[u0,u0+d0*rand(2*N)];diffeq=diffeq)
end

function replace_infs!(Λ,value)
	for i in 1:size(Λ)[1]
		for j in 1:size(Λ)[2]
			if isinf(Λ[i,j])
				Λ[i,j] = value
			end
		end
	end
end

#calculate lyapunov exponents for the phase diagram
#lyapunov considers u0=initial_state(ds)
#system is multistable! final state of p_i should be the initial state of p_{i+1}
#highest level series version
function lyapunov_diagram(ds, gvalues, wvalues; pidx=[1, 2], T, Ttr, Δt, fp_settings=[100, 100, 0.1], d0=1e-9, ut=1e-6, lt=1e-12, parallelize=false, save, realtime_plotting=true, savedir="./")
	gl = length(gvalues)
	wl = length(wvalues)
	Λ = Array{Float64}(undef, (gl,wl))
	fill!(Λ,Inf)
	N = Int(dimension(ds) / 2)
	wc2 = (N-1)/N

	if parallelize
		lyapunov_diagram_parallel!(Λ,ds,gvalues,wvalues; T=T, Ttr=Ttr, Δt=Δt, d0=d0, ut=ut, lt=lt)

	else
		#create pds (parallel system) 
		#for a lower level call to lyapunov
		s = current_state(ds)
		states = [s, s .+ rand(2N) .* d0 / sqrt(2N)]
		pds = ParallelDynamicalSystem(ds,states)

		lyapunov_diagram!(Λ,pds,gvalues,wvalues;pidx=pidx, T=T, Ttr=Ttr, Δt=Δt,fp_settings=fp_settings, d0=d0, ut=ut, lt=lt)
	end

	if save
		replace_infs!(Λ,-1)
		
		gstart = gvalues[1]
		gstop = gvalues[end]
		
		wstart = wvalues[1]
		wstop = wvalues[end]
		
		writedlm(savedir*"N$N"*"_gvalues_gstart_$gstart"*"_gstop_$gstop"*"_T$T"*"_Ttr$Ttr"*"_res_"*"$gl"*"x$wl"*".txt",gvalues)
		writedlm(savedir*"N$N"*"_wvalues_wstart_$wstart"*"_wstop_$wstop"*"_T$T"*"_Ttr$Ttr"*"_res_"*"$gl"*"x$wl"*".txt",wvalues)
		writedlm(savedir*"N$N"*"_T$T"*"_Ttr$Ttr"*"lyaps"*"gvalues_gstart_$gstart"*"_gstop_$gstop"*"_wstart_$wstart"*"_wstop_$wstop"*"_res_"*"$gl"*"x$wl"*".txt",Λ)
		
		if realtime_plotting
			pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...)					
			savefig(pd,savedir*"N$N"*"_T=$T"*"_Ttr=$Ttr"*"_res_"*"$gl"*"x$wl"*"lyap_diagram.pdf")
		end
		
	end 
	return Λ 
end

#parallel lyapunov calculation?
#parallelize row-wise
#each thread saves λ lyap exponent to the same Λ
#pds cannot be shared among threads
function lyapunov_diagram_parallel!(Λ,ds,gvalues,wvalues; T, Ttr, Δt, d0=1e-9, ut=1e-6, lt=1e-12)

	gl = length(gvalues)
	N = Int(dimension(ds)/2)

	#init pds states to initial state of ds
	s = initial_state(ds)
    states = [s, s .+ rand(2N) .* d0 / sqrt(2N)]
	gs_ti = zeros(Int,length(gvalues))

	#start threads wrt g values (columns in the lyap. diag.)
	Threads.@threads for g_i in 1:length(gvalues)
		gs_ti[g_i] = Threads.threadid()
		@show gs_ti
		#@show Threads.threadid(),g_i,gvalues[g_i]
		#create pds (parallel system) for each thread 
		pds = ParallelDynamicalSystem(deepcopy(ds),deepcopy(states))

		set_parameter!(pds,1,gvalues[g_i])
		
		#set direction for w param (lower to larger or viceversa)
		if iseven(g_i)
			witer = reverse(wvalues)
		else
			witer = wvalues
		end
		
		#calc lyapunovs and save to Λ
		#calc for all w values 
		calc_lyap_row!(Λ,pds,gvalues,witer,g_i; T, Ttr, Δt, d0=1e-9, ut=1e-6, lt=1e-12)
	end
end

#serial version (no threading)
function lyapunov_diagram!(Λ,pds,gvalues,wvalues;pidx=[1,2], T, Ttr, Δt,fp_settings=[100,100,0.1] ,d0=1e-9, ut=1e-6, lt=1e-12)

	gl = length(gvalues)

	for (gindex,g) in enumerate(gvalues)
		set_parameter!(pds,pidx[1],g)
		@show g		
		if iseven(gindex)
			witer = reverse(wvalues)
		else
			witer = wvalues
		end
		
		calc_lyap_row!(Λ,pds,gvalues,witer,gindex;pidx=pidx[2], T, Ttr, Δt,fp_settings=fp_settings, d0=d0, ut=ut, lt=lt)
	end
end


function calc_lyap_row!(Λ,pds,gvalues,witer,gindex;pidx=2, T, Ttr, Δt,fp_settings=[100,100,0.1], d0=1e-9, ut=1e-6, lt=1e-12)
	wl = length(witer)
	N = Int(dimension(pds) / 2)
	wc2 = (N-1)/N
	wc1 = (N-1)/gvalues[gindex]
	T_fp,Ttr_fp,Δt_fp = fp_settings

	#reinit!(pds)
	for (windex,w) in enumerate(witer)
		
		#=
		#skip second fixed point parameter interval
		(w <= 1.0 && w > wc2) && continue
		=#

		set_parameter!(pds,pidx,w)
        @show current_parameters(pds)
		reinit!(pds)
		@show w
		
		#parallel debug
		#@show current_state(pds)[1]
		#@show Threads.threadid(),gindex,gvalues[gindex],w

		#use lower transients for stable fixed point dynamics
		if w < wc1 || w > wc2
			lyap = lyapunov(pds, T_fp; Ttr=Ttr_fp, Δt=Δt_fp, d0=d0, d0_upper=ut, d0_lower=lt)
		else
			lyap = lyapunov(pds, T; Ttr=Ttr, Δt=Δt, d0=d0, d0_upper=ut, d0_lower=lt)
		end
		@show lyap
		if iseven(gindex)
			Λ[gindex,wl-windex+1] = lyap
		else
			Λ[gindex,windex] = lyap
		end
	end

end



