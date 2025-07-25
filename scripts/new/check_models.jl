using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics

#convert original params to rescaled
rescaled_params(rp,rc,K,c,α) = ( 1/(1 + rc/(α*K*c*(N-1))), rp/α )

#original params
rp = 1.0;rc = 1.0;K = 1.0;c = 4.0;α = 5*1e-2
p_orig = [rp,rc,K,c,α]

#plot param conversion function
rcs = [1:1:100;];rps = [0.0:0.1:1;];cs = [1:0.1:10;];αs = [1e-2:1e-3:1e-1;]
pl1 = plot(rcs,first.(rescaled_params.(rp,rcs,K,c,α)),xlabel=L"r_c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl2 = plot(rps,last.(rescaled_params.(rps,rc,K,c,α)),xlabel=L"r_p",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl3 = plot(cs,first.(rescaled_params.(rp,rc,K,cs,α)),xlabel=L"c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl4 = plot(αs,first.(rescaled_params.(rp,rc,K,cs,α)),xlabel=L"\alpha",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl5 = plot(αs,last.(rescaled_params.(rp,rc,K,c,αs)),xlabel=L"\alpha",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl_param_conversion = plot(pl1,pl3,pl4,pl5,pl2,layout=(5,1),size=(500,1500),left_margin=8Plots.mm)

#rescaled params
w,g = rescaled_params(rp,rc,K,c,α)
p = [g,w,10]


#initialize systems
u0 = rand(2N)
diffeq = (alg = Rodas5(), reltol = 1e-8, abstol = 1e-10)
ds = create_ds(u0,p;diffeq=diffeq)
ds_orig = create_original_ds(u0,p_orig;diffeq=diffeq)

@time traj,t = trajectory(ds,300;Ttr=2000,Δt=0.01)
@time traj_orig,t_orig = trajectory(ds_orig,5300;Ttr=20000,Δt=0.01)

#timeseries
pl = plot(t,traj[:,1],lc=:gray10,xlabel=L"x_1",ylabel=L"y_1",guidefontsize=20)
pl_orig = plot(t_orig,traj_orig[:,1],lc=:gray10,xlabel=L"L_1",ylabel=L"Y_1",guidefontsize=20)
pl = plot(pl,pl_orig;layout=(2,1),size=(800,400),left_margin=5Plots.mm)

#phase space
pl = plot(traj[:,1],traj[:,2],lc=:gray10,xlabel=L"x_1",ylabel=L"y_1",guidefontsize=20)
pl_orig = plot(traj_orig[:,1],traj_orig[:,2],lc=:gray10,xlabel=L"L_1",ylabel=L"Y_1",guidefontsize=20)
pl = plot(pl,pl_orig;layout=(2,1),size=(400,800),left_margin=5Plots.mm)