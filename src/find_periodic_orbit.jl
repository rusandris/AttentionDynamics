function find_periodic_orbit(att_din,att_din!,prob;init_cond,paramforguess,Ttr,T,dt,M)
		
		#trajectory calculation
		prob_de = ODEProblem(att_din!, init_cond, (0.,Ttr), paramforguess)
		sol = solve(prob_de, Rodas5(),reltol = 1e-8, abstol = 1e-10, adaptive = false, dt = dt)
		prob_de = ODEProblem(att_din!, sol.u[end], (0,T), paramforguess)
		sol = solve(prob_de, Rodas5(),reltol = 1e-8, abstol = 1e-10, adaptive = false, dt = dt)
		
		#period estimation
		x1 = [sol.u[i][1] for i in 1:length(sol.u)]
		time = range(sol.t[1],sol.t[end],length(x1))		
		period = estimate_period(x1,:autocorrelation,time;L=length(x1)÷2,ϵ=0.2)
		
		@info "Estimated period = $period"
		
		#plot ts
		#pl = plot(sol,idxs=[1,2,3,4,5,6],legend=false)
		#display(pl)
		
		
		probtrap = PeriodicOrbitTrapProblem(M = M, N = length(sol.u[end]), prob_vf = prob, xπ = copy(sol[end]), ϕ = copy(sol[end]), updateSectionEveryStep = 0)
		resize!(probtrap.ϕ, probtrap.N*probtrap.M)
		resize!(probtrap.xπ, probtrap.N*probtrap.M)
		
		og = BK.generateSolution(probtrap, t -> sol(t*period/2/π), period)

		_sol = BK.getPeriodicOrbit(probtrap, og, paramforguess)
		
		probtrap.xπ .= 0
		probtrap.ϕ .= reduce(vcat, [att_din(_sol.u[:,i],paramforguess) for i=1:probtrap.M])
		
		solpo = newton(probtrap, og, NewtonPar(verbose = true))
		return probtrap,og
end


function periodic_orbit_continuation(att_din,att_din!;init_cond,paramforguess,Ttr=500.0,T,dt=1e-3,M=300,optnewton = NewtonPar(tol = 1e-11,maxIter = 100 ,verbose = false),pMin = 0.35, pMax = 0.65, dsmax = 0.05, ds , dsmin = 1e-5,jacobian=jacobian,finaliseSolution= (z, tau, step, contResult; prob = nothing, kwargs...) -> true)

	prob = BifurcationProblem(att_din, init_cond, paramforguess, (@lens _.w), J = jacobian,
	recordFromSolution = (x, p) -> norminf(x),
	plotSolution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))

	probtrap,og = find_periodic_orbit(att_din,att_din!,prob;init_cond=init_cond,paramforguess=paramforguess,Ttr=Ttr,T=T,dt=dt,M=M)
	# sup norm
	norminf(x) = norm(x, Inf)
	
	# continuation parameters
	opts_po_cont = ContinuationPar(pMin = pMin, pMax = pMax, dsmax = dsmax, ds = ds, dsmin = dsmin,
 			newtonOptions = optnewton, maxBisectionSteps = 25 ,tolStability = 1e-2,
 			detectBifurcation = 3, saveSolEveryStep=1)


	po = continuation(probtrap, og, PALC(), opts_po_cont;verbosity = 4,
		recordFromSolution = (x, p) -> begin
			N = paramforguess[3];
		 	xtt=reshape(x[1:end-1],2*N,probtrap.M);
		 	fpval = x_fp(p.p,N,1);
		 	fp = zeros(2*N);
		 	fp[[1,N+1]] .= [fpval]
			return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = x[end],
				ymin = minimum(xtt[N+1,:]),
				ymax = maximum(xtt[N+1,:]),
				avg = mean(xtt[1,:]),
				c = centroid(xtt),
				mindist = min_dist(xtt,fp))
			end,
		finaliseSolution = finaliseSolution,
		normC = norminf)
	
end

x_fp(w,N,M) = 1/(1+w*M/(N-1)-w*N/(N-1))

function min_dist(po,fp)
	dim,nr_states = size(po)
	c = zeros(dim)
	min_dist = Inf
	for i in 1:nr_states
		dist_of_state =  norm(po[:,i]-fp)
		if dist_of_state < min_dist
			min_dist = dist_of_state
		end
	end
	return min_dist
end

function centroid(po)
	dim,nr_states = size(po)
	c = zeros(dim)
	for i in 1:nr_states
		c = c + po[:,i]
	end
	return c/nr_states
end


