using AttentionDynamics
using Setfield
using BifurcationKit
using LinearAlgebra
using DelimitedFiles

output_dir = "/home/rusandris/Documents/att_din_project/AttentionDynamics/plotting/data/"

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
par = (g = 5,w = 0.01,N = 3)
norminf(x) = norm(x, Inf)

u0_orbitguess1 = [0.810765070284594, 0.7733829745894634, 0.6229954171573494, 0.5856302648429531, 0.3348947466616211, 0.7743679638708577]
u0_orbitguess2 = [0.2849241076093, 0.08481930770740131, 0.4901061113231271, 0.4619221556280978, 0.27540205357693637, 0.4713032328251774]

initconds = [[0.0,0.0,0.0,0.0,0.0,0.0],
			ones(2*par.N),
			[1/(1-par.w/2),1/(1-par.w/2),0.0,1/(1-par.w/2),1/(1-par.w/2),0.0],
			[1/(1-par.w),0.0,0.0,1/(1-par.w),0.0,0.0]]
			
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
	writedlm(output_dir*"N3fp$i"*"_contdata.txt",hcat(br_fp.param,br_fp.branch.x))
	
	bifpoints = []
	if i == 2 
		push!(bifpoints,[br_fp.specialpoint[1].type, br_fp.specialpoint[1].param,:StableToUnstable])
		push!(bifpoints,[br_fp.specialpoint[2].type, br_fp.specialpoint[2].param,:UnstableToStable])
	
	elseif i == 3
		push!(bifpoints,[br_fp.specialpoint[1].type, br_fp.specialpoint[1].param,:UnstableToUnstable])
		push!(bifpoints,[br_fp.specialpoint[2].type, br_fp.specialpoint[2].param,:UnstableToUnstable])
	
	elseif i != 0
		push!(bifpoints,[br_fp.specialpoint[1].type, br_fp.specialpoint[1].param,:UnstableToStable])
	end
			
	writedlm(output_dir*"N3fp$i"*"_bifpoints.txt",bifpoints) 
	
	
	if i == 3
		# start periodic orbit continuation
		
		optn_po = NewtonPar(tol = 1e-8,  maxIter = 525)
			
		# continuation parameters
		opts_po_cont = ContinuationPar(dsmax = 0.02, ds= 0.001, dsmin = 1e-5,
		 newtonOptions = (@set optn_po.tol = 1e-8), tolStability = 1e-2,
		 detectBifurcation = 3, saveSolEveryStep=1)

		Mt = 500 # number of time sections
		unstable_hopf_po = continuation(
			br_fp, 1, opts_po_cont,
			PeriodicOrbitTrapProblem(M = Mt;
				# specific linear solver for ODEs
				jacobian = :Dense);
			# help branching from Hopf
			usedeflation =true,
			recordFromSolution = (x, p) -> begin
				xtt=reshape(x[1:end-1],6,Mt);
				return (max = maximum(xtt[1,:]),
					min = minimum(xtt[1,:]),
					period = x[end])
				end,
			finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
				# limit the period
				z.u[end] < 1000
				end,
			normC = norminf)
			
		writedlm(output_dir*"N3po_hopf"*"_contdata.txt",hcat(unstable_hopf_po.param,unstable_hopf_po.max))
	end
	
end

po1 = periodic_orbit_continuation(att_din,att_din!;init_cond=u0_orbitguess1,T = 30.0,paramforguess=(g = 5.0,w = 0.41,n = 3),pMin = 0.35,pMax = 0.63,ds = -0.001)

po1_bifpoint_data = []

push!(po1_bifpoint_data,[po1.specialpoint[1].type, po1.specialpoint[1].param,:StableToUnstable]) 
push!(po1_bifpoint_data,[po1.specialpoint[2].type, po1.specialpoint[2].param,:UnstableToStable]) 

writedlm(output_dir*"N3po1_contdata.txt",hcat(po1.param,po1.max))
writedlm(output_dir*"N3po1_contdata_avg.txt",hcat(po1.param,po1.avg))
writedlm(output_dir*"N3po1_bifpoints.txt",po1_bifpoint_data)

po2 = periodic_orbit_continuation(att_din,att_din!;init_cond=u0_orbitguess2,T = 30.0,paramforguess=(g = 5.0,w = 0.41,n = 3),pMin = 0.35,pMax = 0.63,ds = -0.001)
#po1_double = POcont(att_din,att_din!;init_cond=u0_orbitguess1,T = 200.0,paramforguess=(g=5.0,w=0.62,n=3),pMin = 0.61, pMax = 0.65,ds=-0.001)
#po2_double = POcont(att_din,att_din!;init_cond=u0_orbitguess2,T = 200.0,paramforguess=(g=5.0,w=0.62,n=3),pMin = 0.61, pMax =0.64	,ds=-0.001)

po2_bifpoint_data = []

push!(po2_bifpoint_data,[po2.specialpoint[1].type, po2.specialpoint[1].param,:StableToUnstable]) 
push!(po2_bifpoint_data,[po2.specialpoint[2].type, po2.specialpoint[2].param,:UnstableToStable]) 


writedlm(output_dir*"N3po2_contdata.txt",hcat(po2.param,po2.max))
writedlm(output_dir*"N3po2_contdata_avg.txt",hcat(po2.param,po2.avg))
writedlm(output_dir*"N3po2_bifpoints.txt",po2_bifpoint_data)

	
