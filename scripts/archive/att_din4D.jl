using AttentionDynamics
using Setfield
using BifurcationKit
using LinearAlgebra
using DelimitedFiles

function skipsum(a,i,N)
	x = @view a[1:N]
    return sum(x) - x[i]
end


function att_din!(du,u,p,t)
	g,w,n = p
    n = Int(n)
    for i in 1:n
        du[i] = g*u[i]*(1 - (1 - w)*u[i+n] - w/(n-1)*skipsum(u,i,n))
        du[i+n] = u[i] - u[i+n]
    end
    du
end

att_din(u,p) = att_din!(similar(u),u,p,0)

par = (g = 5,w = 0.01,n = 2)
norminf(x) = norm(x, Inf)

initconds = [zeros(2*par.n),ones(2*par.n),[1/(1-par.w),0.0,1/(1-par.w),0.0],[0.0,1/(1-par.w),0.0,1/(1-par.w)]]
#line_colors = [:gray10,:red,:blue,:green]
#labels = [L"$x^{*}_0$",L"$x^{*}_1$",L"$x^{*}_2$",L"$x^{*}_3$"]

for (i,ic) in enumerate(initconds)
	@show ic

	prob = BifurcationProblem(att_din, ic, par, (@lens _.w), J = jacobian,
				recordFromSolution = (x, p) -> norminf(x),
				plotSolution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))


	# continuation options
	opts_br = ContinuationPar(pMin = 0.01, pMax = 0.95, dsmax = 0.01,
		# options to detect bifurcations
		detectBifurcation = 3, nInversion = 8, maxBisectionSteps = 25,
		# number of eigenvalues
		nev = 4, 
		# maximum number of continuation steps
		maxSteps = 1000,
		# parameter theta, see `? continuation`. Setting this to a non
		# default value helps passing the transcritical bifurcation
		θ = 0.3)
	
	br_fp = continuation(prob, PALC(), opts_br;plot = false,normC = norminf)
	writedlm("fp$i"*"_contdata.txt",hcat(br_fp.param,br_fp.branch.x))
	
	bifpoints = []
	if i == 2
		push!(bifpoints,[br_fp.specialpoint[1].type, br_fp.specialpoint[1].param,:StableToUnstable])
		push!(bifpoints,[br_fp.specialpoint[2].type, br_fp.specialpoint[2].param,:UnstableToStable])
	elseif i != 0
		push!(bifpoints,[br_fp.specialpoint[1].type, br_fp.specialpoint[1].param,:UnstableToStable])
	end
			
	writedlm("fp$i"*"_bifpoints.txt",bifpoints) 
	
	
	if i == 2
		# start periodic orbit continuation
		
		optn_po = NewtonPar(tol = 1e-8,  maxIter = 525)
			
		# continuation parameters
		opts_po_cont = ContinuationPar(dsmax = 0.02, ds= 0.001, dsmin = 1e-5,
		 newtonOptions = (@set optn_po.tol = 1e-8), tolStability = 1e-2,
		 detectBifurcation = 3, saveSolEveryStep=1)
		Mt = 500 # number of time sections
		
		br_po = continuation(
			br_fp, 1, opts_po_cont,
			PeriodicOrbitTrapProblem(M = Mt;
				# specific linear solver for ODEs
				jacobian = :Dense);
			# help branching from Hopf
			usedeflation = false,
			recordFromSolution = (x, p) -> (xtt=reshape(x[1:end-1],4,Mt);
				return (max = maximum(xtt[1,:]),
					min = minimum(xtt[1,:]),
					period = x[end])),
			finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
				# limit the period
				z.u[end] < 300
				end,
			normC = norminf)
			
		writedlm("po_hopfN2"*"_bifdata.txt",hcat(br_po.param,br_po.max))
	end
end


