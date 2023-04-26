using BifurcationKit
using Plots,LaTeXStrings
using Setfield
using LinearAlgebra
const BK=BifurcationKit
using DifferentialEquations

include("jacobian.jl")
include("find_periodic_orbit.jl")

norminf(x) = norm(x, Inf)

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

u0_orbitguess = [0.01758938834581525, 0.2604179952706319, 0.4981765183030946, 0.4021602694364159, 0.32005995194379977, 0.4587116869276979, 0.03315556021128763, 0.002418290211501195]

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


w_startvalues_ns = []
w_stopvalues_ns = []
ns_gvalues = [5.5:0.01:9.0;] #g=10 wtf
N=4


for g in ns_gvalues
	wc1 = (N-1)/g
	wc2 = (N-1)/N
	par = (g = g,w = (wc1+wc2)/2,N = N)
	@show par
	po1 = POcont(att_din,att_din!;init_cond=u0_orbitguess,T = 60.0,paramforguess=par,pMin = wc1+0.05, pMax = 0.72,ds = 0.001,finaliseSolution=find_nspoints)
	
	for bp in po1.specialpoint
		@show bp.type
	end
	
	if po1.specialpoint[1].type == :ns && po1.specialpoint[2].type == :ns
		push!(w_startvalues_ns,po1.specialpoint[1].param)
		push!(w_stopvalues_ns,po1.specialpoint[2].param)
	end

end
