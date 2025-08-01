using Plots,LaTeXStrings
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using Attractors
cd(@__DIR__)
figs_dir = "../../paper_figs/"

#traj
diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
w = 0.4
p = [g,w,N]
u0 = rand(2N)
ds = create_ds(u0,p;diffeq=diffeq)

T_limit_cycle = 10
T_plot_trans = 100
Ttr_limc = 1000
Δt = 0.001

#check if two distinct limit cycles can be visualized
#---------------------find two attractors to show as examples----------------------------

cons_recurrences = 500 #nr of recurrences 
locate_steps = 10000 #nr of steps after detection
max_ensemble = 50 #nr of initial conditions

#test trajectory for efficient phase space grid choice
traj,_ = trajectory(ds,T_total;Ttr=500,Δt=1e-3)
Plots.plot(traj[:,1],traj[:,2])
#need max finding along each dim
max_values = maximum(Matrix(traj);dims=1) #ugly max finding
grid_size = 100 #phase space grid cell size

grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])

mapper = AttractorsViaRecurrences(ds, grid;Ttr=1000,
    consecutive_recurrences = cons_recurrences, 
    attractor_locate_steps = locate_steps)

u0s_rand = [rand(2N) for _ in 1:10]
attrs = mapper.(u0s_rand)

attrs_to_plot = [1,2]
#----------------------------plot-------------------------  
colors = [:green,:purple]
pl = plot(xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
        guidefontsize=20,legend=:bottomright,dpi=300)

for (i,attr) in enumerate(attrs_to_plot)
    u0 = u0s_rand[attr]

    #plot transient
    traj_tr,t_tr = trajectory(ds,T_plot_trans,u0;Δt = Δt)

    x1 = traj_tr[:,1]
    x2 = traj_tr[:,3]
    y1 = traj_tr[:,2]
    y2 = traj_tr[:,4]

    plot!(x1 .- x2, y1 .- y2,
        lc=colors[i],lw=1,la=0.2,
        label="")

    #plot limc
    traj_limc,t_limc = trajectory(ds,T_limit_cycle,u0;Ttr=Ttr_limc,Δt = Δt)

    x1 = traj_limc[:,1]
    x2 = traj_limc[:,3]
    y1 = traj_limc[:,2]
    y2 = traj_limc[:,4]

    plot!(x1 .- x2, y1 .- y2,
        lc=colors[i],lw=3,la=1.0,dpi=300,
        label="limit cycle $i",
        title = L"N=%$N"*L", g=%$g"*L", w=%$w")


end
savefig(pl,figs_dir*"two_attractor_diff_proj_test_N$N"*"_g$g"*"_w$w"*".png")

#=
attractors = extract_attractors(mapper)
#plot_attractors(attractors,access=[1,3],framestyle=:auto)
atr1 = attractors[1]
atr2 = attractors[3]
x11 = atr1[:,1] 
x21 = atr1[:,3]
y21 = atr1[:,4] 
y11 = atr1[:,2]
x12 = atr2[:,1]
x22 = atr2[:,3] 
y12 = atr2[:,2]
y22 = atr2[:,4]

plot(atr1[:,1],atr1[:,2],
    lc=:gray80,lw=1,la=1.0,
    xlabel=L"x_1",ylabel=L"y_1",
    guidefontsize=20,legend=false)

plot(x11.-x21,y11.-y21,
    lc=:green,lw=1,la=1.0,
    xlabel=L"x_1-x_2",ylabel=L"y_1-y_2",
    guidefontsize=20,legend=false)

plot!(x12.-x22,y12.-y22,
    lw=1,la=1.0,label="2",
    lc=:purple,
    guidefontsize=20,legend=false)
=#