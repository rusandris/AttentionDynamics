using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
using OrdinaryDiffEq
cd(@__DIR__)

#system def 
N = 10
diffeq = (alg=Tsit5(), dt=1e-3, dtmax=1e-2, adaptive=true, maxiters=1e15)
u0 = rand(2N)
g = 15.75
ds = create_ds(u0, [g, 0.5, N]; diffeq=diffeq) #no noise
ws = [0.5,0.65,0.75,0.92]
Ttrs = [0,1000,1000,0]
#--------------------------------plot time series---------------------------
plots = []
for (i,w) in enumerate(ws)
    set_parameter!(ds, 2, w)
    Ttr = Ttrs[i]
    traj, t = trajectory(ds, 100, u0; Ttr=Ttr, Δt=1e-2)

    pl = plot( xlabel=L"t",
        guidefontsize=18, legendfontsize=9, 
        tickfontsize=10,xlims=(Ttr,Ttr+50))
        #title=L"N = %$N"*L", 1/\epsilon = %$(g), w = %$w")

    plot!(t, Matrix(traj[:,1:2:2N]), lw=3, label="")

    if i == 1
        x_fp = 1.0
        hline!([x_fp], ls=:dash, lc=:gray10,lw=2,label="") #label=L"x^*_N = 1")
    elseif i == 4
        x_fp = 1.0/(1-w)
        hline!([x_fp], ls=:dash, lc=:gray10, lw=2, label="") #label=L"x^*_1 = \frac{1}{1 - w}")
    end
    push!(plots,pl)
end

#adjust plots
plot!(plots[1], ylabel=L"x_i")
x_plots_limits = [(Ttrs[1],Ttrs[1]+50),
    (Ttrs[2],Ttrs[2]+10),(Ttrs[3],Ttrs[3]+10),
    (Ttrs[4],Ttrs[4]+50)]  

for i in 1:length(plots)
    plot!(plots[i],xlims=x_plots_limits[i])
end
plot!(plots[1],left_margin=8Plots.mm)
plot!(plots[4], right_margin=2Plots.mm)

pl4 = plot(plots...,layout=(1,length(plots)),size=(1200,300), dpi = 300,
    bottom_margin=8Plots.mm)

figs_dir = "../../paper_figs/"
savefig(pl4,figs_dir*"timeseries_delayed_LV_N$N"*".png")