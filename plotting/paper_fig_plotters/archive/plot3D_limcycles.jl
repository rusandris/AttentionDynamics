using AttentionDynamics
using OrdinaryDiffEq
using Plots
using LaTeXStrings
plotlyjs()

include("defaults.jl")

p = [4,0.52,3]
ws = [0.52,0.659,0.7]

plots3d = []

for (i,w) in enumerate(ws)
	

	u0 = rand(6)
	p[2] = w 
	g,w,N = p
	Ttr = 100
	Δt = 1e-2
	T = 150
	logocolors = Colors.JULIA_LOGO_COLORS


	ds = create_ds(u0,p)
	diffeq = (alg=Tsit5(),adaptive=false,dt=1e-3)

	state1 = [0.810765070284594, 0.7733829745894634, 0.6229954171573494, 0.5856302648429531, 0.3348947466616211, 0.7743679638708577]
	state2 = [0.2849241076093, 0.08481930770740131, 0.4901061113231271, 0.4619221556280978, 0.27540205357693637, 0.4713032328251774]

	cycle1 = evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = state1,p=p,diffeq = diffeq)
	cycle2= evolve_system(ds,T;Ttr=Ttr,Δt=Δt,state = state2,p=p,diffeq = diffeq)

	pp3d = plot()


	x3(x1,x2) = 3 - x1 - x2
	x1 = [0:0.01:2.0;]
	x2 = [0:0.01:2.0;]

	surface!(x1,x2,x3;THREED_DEFAULTS...,camera=(-30,30),color=:red,alpha=0.4,colorbar=false,size = (1000, 1000),grid=true)

	x_fp(w,N,M) = 1/(1+w*M/(N-1)-w*N/(N-1))
	fp_c = [1.0,1.0,1.0]
	fp_x1 = [x_fp(w,3,1),0,0]
	fp_x2 = [0,x_fp(w,3,1),0]
	fp_x3 = [0,0,x_fp(w,3,1)]

	plot3d!([fp_c[1],fp_x1[1],fp_x2[1],fp_x3[1]],[fp_c[2],fp_x1[2],fp_x2[2],fp_x3[2]],[fp_c[3],fp_x1[3],fp_x2[3],fp_x3[3]],st=:scatter,mc=:gray10,ma=1,ms=4,label="fp")
	
	if i != 3		
		plot3d!(pp3d,cycle1[:,1],cycle1[:,3],cycle1[:,5];THREED_DEFAULTS...,lc=:gray10,xlabel="x1",ylabel="x2",zlabel="x3")
		plot3d!(pp3d,cycle2[:,1],cycle2[:,3],cycle2[:,5];THREED_DEFAULTS...,lc=:gray50,xlabel="x1",ylabel="x2",zlabel="x3")
	end
		
	push!(plots3d,pp3d)



end

plot(plots3d...,layout=(1,3),size=(1000,500))

