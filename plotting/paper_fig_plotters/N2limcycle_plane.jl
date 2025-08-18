using AttentionDynamics
using OrdinaryDiffEq
using Plots
using LaTeXStrings
using LinearAlgebra
using DynamicalSystemsBase
using DelimitedFiles
cd(@__DIR__)
include("../defaults.jl")
include("../plot_functions/timeseries_phase_portrait.jl")

#include("bif_floq_od_method_parameters.jl")
#run necessary mkpaths
#create dir for figs
figs_dir = "../../paper_figs/"
mkpath(figs_dir)

#----------define local parameters------------
N=2 #small system
Ttr_N2 = 0
Δt = 1e-3


#--------------------------plotting params-------------------------
logocolors = Colors.JULIA_LOGO_COLORS
w_values = [0.3,0.41]
Ttrs = [200,200]
Ts = [20,200]
pp_lims = [(0.0,3.0),(0.0,4.0)]
pp_ticks = [[0.0,1.5,3.0],[0.0,2.0,4.0]]
ts_ticks = [[0.0,1.5,3.0],[0.0,2.5,5.0]]
ts_lims = [(0.0,3.0),(0.0,5.0)]
velm_lims = [(1.0,4.0),(0.0,5.5)]
velm_ticks = [[1.0,2.5,4.0],[0.0,2.5,5.5]]


diffeq = (alg = Rodas5(),dt = 1e-3,adaptive=true,reltol = 1e-8, abstol = 1e-10,maxiters=1e15)

g = 5.0
wc1 = (N-1)/g
wc2 = (N-1)/N

#generate random initial conditions
u0s = [rand(2N) for i in 1:3]
#u0 = generate_initial_condition(N)

p = [g,w_values[1],N]
@show wc1,wc2,p
#create system
ds = create_ds(u0s[1],p;diffeq=diffeq)
colors = [logocolors.blue, logocolors.green, :gray50]
plots = []

#plot for two w values 
for i in 1:length(w_values)
	w = w_values[i]
    @show w
	lims=pp_lims[i]
	ticks=round.(pp_ticks[i];digits=1)
	velm_yticks = velm_ticks[i]
	ts_ylims=ts_lims[i]
	velm_ylims=velm_lims[i]
	ts_yticks = ts_ticks[i] 
	set_parameter!(ds,2,w)
	
    pp = plot(xlabel=L"$x_1$",ylabel=L"$x_2$",
        lims=lims,la=0.9,ticks=ticks,
        title=L"w=%$w",titlefontsize=15)

    if i == 1
        plot!(xformatter=:none,xlabel="")
    end

	#calc trajectories for all init cond
    for (j,u0) in enumerate(u0s)
        traj,t = trajectory(ds,Ttrs[i],u0;Ttr=0,Δt = 1e-3)
    
        #plot phase portrait projection x1, x2
        plot!(pp,traj[:,1],traj[:,3];DEFAULTS...,
            lw=2,la=0.4,lc=colors[j],label="")

        #plot stable cycle as well
        cycle,t = trajectory(ds,Ts[i],traj[end];Ttr=Ttr_N2,Δt = 1e-3)

        #plot phase portrait projection x1, x2
        plot!(pp,cycle[:,1],cycle[:,3];DEFAULTS...,
            lw=4,la=1.0,lc=colors[j],label="")

    end

	#---------------degeneracy plane------------------
	x2(x1) = 2 - x1 
	plot!(pp,[0,2],x2.([0,2]),lc=:gray10,la=0.8,ls=:dash,lw=5,label="degeneracy plane")

	#---------------------fixed points----------------------
	#dominance fixed points
	u_d1 = zeros(2N)
	u_d2 = zeros(2N)

	#calc dominance fp value
	M=1
	w_M = (w*M)/(N-1)
	w_N = (w*N)/(N-1)
	s = 1/(1 + w_M - w_N) #fixed point value (nonzero)

	u_d1[1:2] .= [s,s] #set both x1 and y1 to s
	u_d2[3:4] .= [s,s] #set both x2 and y2 to s

	@show u_d1,u_d2
	# coexistence fp
	u_c = ones(2N)
	plot!(pp,[u_c[1]],[u_c[2]],st=:scatter,mc=:white,markerstrokecolor=:red,markerstrokewidth=2,ma=1,ms=8,label=L"s_c^*",legend=true,legendfontsize=12)

	#dominance fp
	plot!(pp,[u_d1[1],u_d2[1]],[u_d1[3],u_d2[3]],st=:scatter,mc=:white,markerstrokecolor=:blue,markerstrokewidth=2,ma=1,ms=8,label=L"s_d^*",legend=true,legendfontsize=12)

    push!(plots,pp)

	#combine all plots
	#pl = plot(pp,ts_comb,layout=(1,2),size=fig_size,bottom_margin=10Plots.mm,left_margin=10Plots.mm)
	#savefig(pl,figs_dir*"timeseries_phase_portrait_N2_degen_plane_w=$w"*".pdf")

end
pl_limc = plot(plots...,layout=(2,1),size=(500,800),bottom_margin=3Plots.mm,left_margin=3Plots.mm)
savefig(pl_limc,figs_dir*"phase_portrait_N2_g$(g)_degen_plane"*".pdf")


#=
#----------------------phase space velocity---------------------------
velocity = zeros(4)
velocity_magnitudes = zeros(length(cycle[:,1]))

for i in 1:size(cycle)[1]
    att_din!(velocity,cycle[i],p,0) #in-place call to calculate derivatives for every state
    velocity_magnitudes[i] = norm(velocity)
end

#plot ts
ts = plot(t[:],Matrix(cycle[:,[1,3]]);lc=[:gray10 :gray50],lw=3,
    ms=0,mc=:purple,markerstrokewidth=0.0,
    ls=:solid,xformatter=:none, ylabel=L"x_i",
    labels=[L"x_1" L"x_2"],DEFAULTS...)

#plot velocity 
ts_vel = plot(t[:],velocity_magnitudes[:];lc=:gray10,lw=3,
    markerstrokewidth=0.0,legend=false,
    ylabel=L"$\left|\dot{\mathbf{s}}\right|$",xlabel=L"t",DEFAULTS...)

#combine the two
ts_comb = plot(ts,ts_vel,layout=(2,1))
=#