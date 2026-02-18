using DynamicalSystems,AttentionDynamics,OrdinaryDiffEq,DelimitedFiles,Plots,LaTeXStrings

include("defaults.jl")
include("plot_functions/timeseries_phase_portrait.jl")


ics = readdlm("data/N4_LIMCYCLES_IC.txt")

u0 = ics[:,5]
special_params = [(w = 0.7,g = 6.0),(w = 0.71,g = 6.0),(w = 0.692,g = 8.0)]

p = [6.0,0.67,4]
ds = create_ds(u0,p)

T = 500
Ttr = 3000
Δt = 1e-2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=false,maxiters=1e15)
N = 4

pp_colors = [:gray70, :gray40, :gray10]
pp_labels = (x = L"x_1",y =  L"y_1")
ts_labels = (x = L"t",y = ylabel = L"$x_i$")

#pp_labels = (x = "",y =  "")
#ts_labels = (x = "",y ="")

all_plots = []

for (i,params) in enumerate(special_params)
	g = params.g
	w = params.w
	@show g,w
	set_parameter!(ds,1,g)
	set_parameter!(ds,2,w)
	
	ts = evolve_system(ds,T;Ttr = Ttr,Δt = 1e-2,state = u0,p,diffeq=diffeq)
	xts = ts[:,1]
	yts = ts[:,2]

	if i == 3 
		ts_xlabel = ts_labels.x
		pp_xlabel = pp_labels.x
	else
		ts_xlabel = ""
		pp_xlabel = ""
	end
	
	ts_plot_height = 12
	pp_plot_height = 8
	pp_plot_ticksmax = maximum(ts[:,1])+1

	plot_ts = plot_timeseries(ts;DEFAULTS...,T=T,Ttr=Ttr,Δt=Δt,lw=1.5,
	xlims = (Ttr,Ttr+30),ylims=(0,ts_plot_height),widen=true,
	ylabel=ts_labels.y,xlabel=ts_xlabel,yticks=0:5:ts_plot_height,right_margin=10Plots.mm)
	
	#display(plot_ts)
	
	plot_pp = plot_phase_portrait(ts[:,1],ts[:,2];DEFAULTS...,lc=pp_colors[i],legend=true,label="w = $w"*" g = $g",xlims = (0,maximum(xts)),ylims = (0,pp_plot_height),yticks = 0:4:pp_plot_height,xticks = 0:4:pp_plot_ticksmax,widen=true,xlabel=pp_xlabel,ylabel=pp_labels.y)
	
	push!(all_plots,plot_ts)
	push!(all_plots,plot_pp)
	
	#savefig(plot_ts,"ts_N$N"*"_w$w"*"_g$g"*"_T$T"*"_Ttr$Ttr"*".png")
	#savefig(plot_pp,"pp_N$N"*"_w$w"*"_g$g"*"_T$T"*"_Ttr$Ttr"*".png")
end
l = @layout [a{0.6w} b{0.2w}; c{0.6w} d{0.2w}; e{0.6w} f{0.2w}]
all = plot(all_plots...,layout=l)

savefig(all,"figs/all_ts_pp.pdf")
