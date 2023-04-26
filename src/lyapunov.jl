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



function lyapunov_diagram(pinteg,gvalues,wvalues; T, Ttr, Δt, d0=1e-9, ut=1e-6, lt=1e-12,save,realtime_plotting=true)
	gl = length(gvalues)
	wl = length(wvalues)
	Λ = Array{Float64}(undef, (gl,wl))
	fill!(Λ,Inf)
	N = pinteg.p[3]
	wc2 = (N-1)/N

	for (gindex,g) in enumerate(gvalues)
		set_parameter!(pinteg,1,g)
		
		if iseven(gindex)
			witer = reverse(wvalues)
		else
			witer = wvalues
		end
		
		for (windex,w) in enumerate(witer)
			
			#skip second fixed point parameter interval
			(w <= 1.0 && w > wc2) && continue
	
			set_parameter!(pinteg,2,w)
			reinit!(pinteg)

			@show (g,w)
			lyap = lyapunov(pinteg, T, Ttr, Δt, d0, ut, lt)
			@show lyap
			if iseven(gindex)
				Λ[gindex,wl-windex+1] = lyap
			else
				Λ[gindex,windex] = lyap
			end
			
			if realtime_plotting
				pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...)					
				display(pd)
			end
			
		end
	end
	
	if save
		replace_infs!(Λ,-1)
		
		writedlm("N=$N"*"_gvalues"*"_T=$T"*"_Ttr=$Ttr"*"_res_"*"$gl"*"x$wl"*".txt",gvalues)
		writedlm("N=$N"*"_wvalues"*"_T=$T"*"_Ttr=$Ttr"*"_res_"*"$gl"*"x$wl"*".txt",wvalues)
		writedlm("N=$N"*"_T=$T"*"_Ttr=$Ttr"*"lyaps"*"_res_"*"$gl"*"x$wl"*".txt",Λ)
		
		if realtime_plotting
			pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...)					
			savefig(pd,"N=$N"*"_T=$T"*"_Ttr=$Ttr"*"_res_"*"$gl"*"x$wl"*"lyap_diagram.pdf")
		end
		
	else 
		return Λ 
	end
end

