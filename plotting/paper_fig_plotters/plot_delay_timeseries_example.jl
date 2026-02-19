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

#---------------------read parameter values------------------
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c/"
params_orig = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled = readdlm(param_dir * "params_w_g.txt")
params_changed = readdlm(param_dir * "params_c.txt")
N_orig = params_orig[end]

#=
w_B2 = 0.48
ws = [0.3, w_B2, 0.7, 0.95]
g = 22.5 #15.75
=#
ws = params_rescaled[:, 1]
gs = params_rescaled[:, 2] ./ (N_orig - 1) * ((N - 1))
g = gs[1]

ds = create_ds(u0, [g, 0.5, N]; diffeq=diffeq) #no noise

#ws = [0.5,0.65,0.75,0.92]
Ttrs = [0,1000,1000,0]
T = 200

#--------------------------------plot time series---------------------------
#check WTA 
set_parameter!(ds,2,ws[end])
traj_WTA, t = trajectory(ds, T, u0; Ttr=Ttrs[end], Δt=1e-2)
x_idx_WTA, = findall(x->x>1.0,traj_WTA[end])

plots = []
colors = collect(Colors.JULIA_LOGO_COLORS)[[1, 2, 4]]  #use three colors
x_idxs = [1, 3, x_idx_WTA]


for (i,w) in enumerate(ws)
    set_parameter!(ds, 2, w)
    Ttr = Ttrs[i]
    traj, t = trajectory(ds, T, u0; Ttr=Ttr, Δt=1e-2)

    pl = plot( xlabel=L"t",
        guidefontsize=18, legendfontsize=9, 
        tickfontsize=10,xlims=(Ttr,Ttr+50))
        #title=L"N = %$N"*L", 1/\epsilon = %$(g), w = %$w")

    #plot some with gray
    plot!(t, Matrix(traj[:,7:2:20]), lw=1, lc=:gray10, la=0.3, label="")

    #plot 3 with colors 
    plot!(t, Matrix(traj[:, x_idxs]), lw=3, lc=colors', label="")

    if i == 1
        x_fp = 1.0
        hline!([x_fp], ls=:dash, lc=:gray10,lw=2,label="") #label=L"x^*_N = 1")
    elseif i == 4
        x_fp = 1.0/(1-w)
        hline!([x_fp], ls=:dash, lc=:gray10, lw=2, label="") #label=L"x^*_1 = \frac{1}{1 - w}")
    end
    push!(plots,pl)
end


x_plots_limits = [(Ttrs[1],Ttrs[1]+50),
                (Ttrs[2],Ttrs[2]+20),
                (Ttrs[3],Ttrs[3]+40),
                (Ttrs[4],Ttrs[4]+20)]
                
                
y_plots_limits = [(0,4),(0,12),(0,40),(0,100)]

#adjust time lims and ticks 
xticks = [[0,25,50],
    [1000, 1010, 1020],
    [1000, 1020, 1040],
    [0, 10 ,20]
    ]

yticks = [[0,1,2,3,4],
    [0, 4, 8, 12],
    [0,10,20,30],
    [0, 40 ,80]
    ]

for i in 1:length(plots)
    plot!(plots[i],xlims=x_plots_limits[i],
        xticks=xticks[i], yticks=yticks[i], ylims=y_plots_limits[i])
end

#adjust plots
plot!(plots[1], ylabel=L"x_i",left_margin=5Plots.mm)
plot!(plots[3], ylabel=L"x_i", left_margin=5Plots.mm)
plot!(plots[1],xlabel="")
plot!(plots[2],xlabel="",right_margin=5Plots.mm)

#=
#add annotations 
labels = [L"A",L"B_2",L"C",L"D"]
x_annot_offsets = [-0.5,-0.5,-0.5,-0.5]
for (i,pl) in enumerate(plots)
    xpos = x_plots_limits[i][1]+x_annot_offsets[i]
    ypos = y_plots_limits[i][end]*0.9
    annotate!(pl,xpos,ypos, text(labels[i], :gray10, :center, 18))
end
=#
pl4 = plot(plots...,layout=(2,2),size=(800,500), dpi = 300)

figs_dir = "../../paper_figs/"
savefig(pl4,figs_dir*"timeseries_delayed_LV_N$(N)_g$(g)"*"labeled_threecolors_2x2.png")