using BifurcationKit
using Setfield
using LinearAlgebra
const BK=BifurcationKit
using DifferentialEquations
using StatsBase
using AttentionDynamics
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



u0_orbitguess1 = [0.810765070284594, 0.7733829745894634, 0.6229954171573494, 0.5856302648429531, 0.3348947466616211, 0.7743679638708577]
u0_orbitguess2 = [0.2849241076093, 0.08481930770740131, 0.4901061113231271, 0.4619221556280978, 0.27540205357693637, 0.4713032328251774]

init_conds = [u0_orbitguess1,u0_orbitguess2,u0_orbitguess1,u0_orbitguess2]
og_par = (g = 5,w = 0.41,n = 3)

po1 = periodic_orbit_continuation(att_din,att_din!;init_cond=u0_orbitguess1,T = 30.0,paramforguess=(g = 5.0,w = 0.41,n = 3),pMin = 0.35,pMax = 0.63,ds = -0.001)
po2 = periodic_orbit_continuation(att_din,att_din!;init_cond=u0_orbitguess2,T = 30.0,paramforguess=(g = 5.0,w = 0.41,n = 3),pMin = 0.35,pMax = 0.63,ds = -0.001)
#po1_double = POcont(att_din,att_din!;init_cond=u0_orbitguess1,T = 200.0,paramforguess=(g=5.0,w=0.62,n=3),pMin = 0.61, pMax = 0.65,ds=-0.001)
#po2_double = POcont(att_din,att_din!;init_cond=u0_orbitguess2,T = 200.0,paramforguess=(g=5.0,w=0.62,n=3),pMin = 0.61, pMax =0.64	,ds=-0.001)

writedlm("po1_bifdata.txt",hcat(po1.param,po1.mindist))
writedlm("po2_bifdata.txt",hcat(po2.param,po2.mindist))

bifpoint_data = []

push!(bifpoint_data,[po1.specialpoint[1].type, po1.specialpoint[1].param,-1]) 
push!(bifpoint_data,[po1.specialpoint[2].type, po1.specialpoint[2].param,1]) 

writedlm("po1_bifpoints.txt",bifpoint_data)

