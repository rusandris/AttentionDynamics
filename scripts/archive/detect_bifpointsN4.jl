using AttentionDynamics
using LinearAlgebra
using DelimitedFiles
norminf(x) = norm(x, Inf)

#-------------------------define system-----------------------------
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
par = (g = 6,w = 0.1,N = 4)

u0s = readdlm("N4_LIMCYCLES_IC.txt")
ic_index = 5
u0_orbitguess = u0s[:,ic_index]

#=
po1 = POcont(att_din,att_din!;init_cond=u0_orbitguess,T = 30.0,paramforguess=(g = 6.0,w = 0.51,n = 4),pMin = 0.45, pMax = 0.73,ds = -0.001)

optnewton = NewtonPar(tol = 1e-11,maxIter = 100 ,verbose = false)
opts_codim2_cont = ContinuationPar(pMin = 3.0, pMax = 10.0, dsmax = 0.05, ds = 0.01, dsmin = 1e-3,
		newtonOptions = optnewton, maxBisectionSteps = 25 ,tolStability = 1e-2,
		detectBifurcation = 3, saveSolEveryStep=1)
=#

function find_nspoints(z, tau, step, contResult;prob = nothing, kwargs...)
	if	length(contResult.specialpoint) == 2
		if contResult.specialpoint[1].type == :ns && contResult.specialpoint[2].type == :ns
			return false
		else
			return true
		end
	else return true
	end
end


w_starts = []
w_stops = []
gvalues = [5.0:0.1:10.0;] #g=10 wtf
N=4


for g in gvalues
	wc1 = (N-1)/g
	wc2 = (N-1)/N
	par = (g = g,w = (wc1+wc2)/2,N = N)
	@show par
	cont_res = periodic_orbit_continuation(att_din,att_din!;init_cond=u0_orbitguess,T = 60.0,paramforguess=par,pMin = wc1+0.05, pMax = 0.72,ds = 0.001)
	
	for bp in cont_res.specialpoint
		@show bp.type
		bptype = bp.type
		wval = bp.param
		
		if bptype != :endpoint 
		
			open("N4bifpoints.txt","a") do f
				write(f,"$g \t"*"\t $bptype"*"\t $wval\n")
			end
		end
		
		#push!(w_starts,cont_res.specialpoint[1].param)
		#push!(w_stops,cont_res.specialpoint[2].param)
	end

end
