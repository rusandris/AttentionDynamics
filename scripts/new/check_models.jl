using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics
using OrdinaryDiffEq

#convert original params to rescaled
rescaled_params(rp,rc,K,c,α,N) = ( 1/(1 + rc/(α*K*c*(N-1))), rp/α )

#original params
N = 10
rp = 1.0;rc = 1.0;K = 1.0;c = 2.0;α = 5*1e-2
p_orig = [rp,rc,K,c,α]

#plot param conversion function
rcs = [1:1:100;];rps = [0.0:0.1:1;];cs = [1:0.1:10;];αs = [1e-2:1e-3:1e-1;]
pl1 = plot(rcs,first.(rescaled_params.(rp,rcs,K,c,α,N)),xlabel=L"r_c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl2 = plot(rps,last.(rescaled_params.(rps,rc,K,c,α,N)),xlabel=L"r_p",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl3 = plot(cs,first.(rescaled_params.(rp,rc,K,cs,α,N)),xlabel=L"c",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl4 = plot(αs,first.(rescaled_params.(rp,rc,K,cs,α,N)),xlabel=L"\alpha",ylabel=L"w",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl5 = plot(αs,last.(rescaled_params.(rp,rc,K,c,αs,N)),xlabel=L"\alpha",ylabel=L"g",
    lw=2,lc=:gray10,legend=false,guidefontsize=20)

pl_param_conversion = plot(pl1,pl3,pl4,pl5,pl2,layout=(5,1),size=(500,1500),dpi=300,left_margin=8Plots.mm)
figs_dir = "../../paper_figs/"
savefig(pl_param_conversion,figs_dir*"pl_param_conversion_test.png")

#rescaled params
w,g = rescaled_params(rp,rc,K,c,α,N)
p = [g,w,10]

#initialize systems
u0 = rand(2N)
diffeq = (alg = Rodas5(), reltol = 1e-8, abstol = 1e-10)
ds = create_ds(u0,p;diffeq=diffeq)
ds_orig = create_original_ds(u0,p_orig;diffeq=diffeq)

@time traj,t = trajectory(ds,300;Ttr=2000,Δt=0.01)
@time traj_orig,t_orig = trajectory(ds_orig,5300;Ttr=20000,Δt=0.01)

#timeseries
pl = plot(t,traj[:,1],lc=:gray10,xlabel=L"t",ylabel=L"x_1",guidefontsize=20)
pl_orig = plot(t_orig,traj_orig[:,1],lc=:gray10,xlabel=L"t",ylabel=L"L_1",guidefontsize=20)
pl_ts = plot(pl,pl_orig;layout=(2,1),size=(800,400),left_margin=5Plots.mm,legend=false,
    plot_title = L"r_p = %$(rp), r_c = %$(rc), K = %$(K),c = %$(c),\alpha = %$(α)"  )

figs_dir = "../../paper_figs/"

savefig(pl_ts,figs_dir*"orig_vs_modelled_ts_c_$c"*".png")

#phase space
pl = plot(traj[:,1],traj[:,2],lc=:gray10,xlabel=L"x_1",ylabel=L"y_1",guidefontsize=20,legend=false)
pl_orig = plot(traj_orig[:,1],traj_orig[:,2],lc=:gray10,xlabel=L"L_1",ylabel=L"Y_1",guidefontsize=20,legend=false)
pl_pp = plot(pl,pl_orig;layout=(2,1),size=(400,800),left_margin=5Plots.mm,
    plot_title = L"r_p = %$(rp), r_c = %$(rc)," * "\n" * L" K = %$(K),c = %$(c),\alpha = %$(α)" )

savefig(pl_pp,figs_dir*"orig_vs_modelled_pp_c_$c"*".png")

