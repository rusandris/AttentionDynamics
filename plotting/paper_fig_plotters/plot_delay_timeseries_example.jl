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
                (Ttrs[2],Ttrs[2]+10),
                (Ttrs[3],Ttrs[3]+20),
                (Ttrs[4],Ttrs[4]+20)]
                
                
y_plots_limits = [(0,6),(0,16),(0,40),(0,100)]

#adjust time lims and ticks 
ticks = [[0,25,50],
    [1000, 1005, 1010],
    [1000, 1010, 1020],
    [0, 10 ,20]
    ]
for i in 1:length(plots)
    plot!(plots[i],xlims=x_plots_limits[i],
        xticks=ticks[i], ylims=y_plots_limits[i])
end
plot!(plots[1],xlabel="")
plot!(plots[2],xlabel="",right_margin=5Plots.mm)

#add annotations 
labels = ["A","B","C","D"]
x_annot_offsets = [5,0.5,1,1]
for (i,pl) in enumerate(plots)
    xpos = x_plots_limits[i][1]+x_annot_offsets[i]
    ypos = y_plots_limits[i][end]*0.9
    annotate!(pl,xpos,ypos, text(labels[i], :gray10, :center, 18))
end

pl4 = plot(plots...,layout=(2,2),size=(600,600), dpi = 300)

figs_dir = "../../paper_figs/"
savefig(pl4,figs_dir*"timeseries_delayed_LV_N$N"*"labeled2x2.png")