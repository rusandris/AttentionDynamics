using StateTransitionNetworks
using DynamicalSystems  
using Plots,LaTeXStrings
cd(@__DIR__)

#chaotic
ds = Systems.logistic()
set_parameter!(ds, 1, 3.82)
traj, t = trajectory(ds, 100; Ttr=1000)

#periodic
set_parameter!(ds,1,3.835)
traj_per, t = trajectory(ds, 100; Ttr=1000)

#random
ts_rand = rand(100)

#plots 

pl1 = plot(t, traj[:, 1], markershape=:circle,
    lc=:gray10, mc=:gray10, ylims=(0, 1.2),
    xlabel=L"\tau", ylabel=L"x_{\tau}", label=L"r=4.0",
    guidefontsize=18, legendfontsize=12)

pl2 = plot(t, traj_per[:, 1], markershape=:circle,
    lc=:gray10, mc=:gray10, ylims=(0, 1.2),
    xlabel=L"\tau", ylabel=L"x_{\tau}", label=L"r=3.835",
    guidefontsize=18, legendfontsize=12)

pl3 = plot(ts_rand, markershape=:circle,
    lc=:gray10, mc=:gray10, ylims=(0, 1.2),
    xlabel=L"\tau", ylabel=L"x_{\tau}", label="random uniform",
    guidefontsize=18, legendfontsize=12,legend=:topright)

#long trajectories
set_parameter!(ds, 1, 3.82)
traj, t = trajectory(ds, 2000; Ttr=1000)
set_parameter!(ds, 1, 3.835)
traj_per, t = trajectory(ds, 2000; Ttr=1000)
ts_rand = rand(2000)

#measures 
sts = timeseries_to_grid(traj, 10; grid_edges=[0.0, 1.0])
P, _, x = calculate_transition_matrix(sts)
@show network_measures(P)

sts_per = timeseries_to_grid(traj_per, 10; grid_edges=[0.0, 1.0])
P_per, _, x = calculate_transition_matrix(sts_per)
@show network_measures(P_per)

sts_rand = timeseries_to_grid(ts_rand, 10; grid_edges=[0.0, 1.0])
P_rand, _, x_rand = calculate_transition_matrix(sts_rand)
@show network_measures(P_rand)



pp1 = plot(traj[1:end-1, 1], traj[2:end, 1], markershape=:circle,
    st=:scatter, mc=:gray10, ms=1.5, ma=0.01,
    ylims=(0, 1.2), xlims=(0, 1.2),
    xlabel=L"x_\tau", ylabel=L"x_{\tau+1}", label="",
    guidefontsize=18, legendfontsize=12)

pp2 = plot(traj_per[1:end-1, 1], traj_per[2:end, 1], markershape=:circle,
    st=:scatter, mc=:gray10, ylims=(0, 1.2), xlims=(0, 1.2), ms=5, ma=1,
    xlabel=L"x_\tau", ylabel=L"x_{\tau+1}", label="",
    guidefontsize=18, legendfontsize=12)



pp3 = plot(ts_rand[1:end-1, 1], ts_rand[2:end, 1], markershape=:circle,
    st=:scatter, mc=:gray10, ms=0.5, ma=0.7,
    ylims=(0, 1.2), xlims=(0, 1.2),
    xlabel=L"x_\tau", ylabel=L"x_{\tau+1}", label="",
    guidefontsize=18, legendfontsize=12)

l = @layout [a{0.6w} a{0.4w}; a{0.6w} a{0.4w}; a{0.6w} a{0.4w}]
pl = plot(pl2, pp2, pl1, pp1, pl3, pp3, layout=l, size=(800, 600),
    left_margin=5Plots.mm, dpi=300)

pl1_resized = plot(pl1,size=(800,300),left_margin=5Plots.mm,bottom_margin=5Plots.mm, dpi=300)
savefig(pl1_resized,"logistic_rand_r4.png")
savefig(pl,"logistic_rand.png")