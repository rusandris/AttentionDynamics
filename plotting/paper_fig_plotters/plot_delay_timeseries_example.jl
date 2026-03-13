using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

#system def 
N = 10
diffeq = (alg=Tsit5(), dt=1e-2, dtmax=1e-3, adaptive=true, maxiters=1e15)
u0 = rand(2N)

#---------------------read parameter values------------------

#horizontal direction on phasediagram 
param_dir_root = "../../data/parameters/"
param_dir = param_dir_root * "c_article/"
params_orig_horizontal = readdlm(param_dir * "params_rp_rc_K_α_N.txt")
params_rescaled_horizontal = readdlm(param_dir * "params_w_g.txt")
params_changed_horizontal = readdlm(param_dir * "params_c.txt")
params_changed_horizontal = vec(params_changed_horizontal)
cs = vec(params_changed_horizontal)


#alpha (rp=rc) direction on phasediagram
param_dir = param_dir_root * "alpha_rp=rc_article/"
params_orig_curve = readdlm(param_dir * "params_K_c_alpha_N.txt")
K, c, α, N_orig = params_orig_curve
params_rescaled_curve = readdlm(param_dir * "params_w_g.txt")
params_changed_curve = readdlm(param_dir * "params_r.txt")
rs = vec(params_changed_curve)

#equal distance to coexistence curve
#on g/(N-1) graph
#between big and small systems 

N_orig = params_orig_horizontal[end]
N_small = 10
ws_horizontal = params_rescaled_horizontal[:, 1]
gs_horizontal = params_rescaled_horizontal[:, 2] * (N_small - 1)
ws_curve = params_rescaled_curve[:, 1]
gs_curve = params_rescaled_curve[:, 2] * (N_small - 1)

#------------------parameters for example timeseries-----------------
g = gs_horizontal[1]
ws_example = [0.1, ws_horizontal[2], ws_horizontal[3], 0.95]

#-------------------define system-------------------

ds = create_ds(u0, [g, 0.5, N]; diffeq=diffeq) #no noise

#ws = [0.5,0.65,0.75,0.92]
Ttrs = [0,100,100,0]
T = 200

#--------------------------------plot time series---------------------------
#check WTA 
set_parameter!(ds,2,ws_example[end])
traj_WTA, t = trajectory(ds, 0.0, u0; Ttr=200, Δt=1e-2)
x_idx_WTA = findfirst(x->x>1.0,traj_WTA[end])

plots = []
colors = collect(Colors.JULIA_LOGO_COLORS)[[1, 2, 4]]  #use three colors
titles = [L"A",L"B_2",L"C",L"D"]
title_colors = [:orange :lime :gray15 :gray80]


guidefontsize = 24
legendfontsize = 12
tickfontsize = 14
titlefontsize = 24

if x_idx_WTA != 1 && x_idx_WTA != 19
    x_idxs = [x_idx_WTA-2, x_idx_WTA+2, x_idx_WTA] #plot with colors
elseif x_idx_WTA == 1
    x_idxs = [x_idx_WTA + 2, x_idx_WTA + 4, x_idx_WTA]
elseif  x_idx_WTA == 19
    x_idxs = [x_idx_WTA - 2, x_idx_WTA - 4, x_idx_WTA]
end

#plot with gray
x_idxs_others = []
for i in 1:2:2N
    i ∈ x_idxs && continue
    push!(x_idxs_others, i)
end

#loop through the parameters
for (i,w) in enumerate(ws_example)
    set_parameter!(ds, 2, w)
    Ttr = Ttrs[i]
    traj, t = trajectory(ds, T, u0; Ttr=Ttr, Δt=1e-2)

    pl = plot( xlabel=L"t",title=titles[i],
        guidefontsize=guidefontsize, legendfontsize=legendfontsize,
        tickfontsize=tickfontsize, titlefontsize=titlefontsize, 
        titlefontcolor=title_colors[i],legend=false, xlims=(Ttr, Ttr + 50))
        #title=L"N = %$N"*L", 1/\epsilon = %$(g), w = %$w")

    #plot some with gray
    plot!(t, Matrix(traj[:,x_idxs_others]), lw=1, lc=:gray10, la=0.3, label="")

    #plot 3 with colors 
    plot!(t, Matrix(traj[:, x_idxs]), lw=2, lc=colors', label="")

    if i == 1
        x_fp = 1.0
        hline!([x_fp], ls=:dash, lc=:gray10,lw=2,label="") #label=L"x^*_N = 1")
    elseif i == 4
        x_fp = 1.0/(1-w)
        hline!([x_fp], ls=:dash, lc=:gray10, lw=2, label="") #label=L"x^*_1 = \frac{1}{1 - w}")
    end
    push!(plots,pl)
end


x_plots_limits = [(Ttrs[1],Ttrs[1]+40),
                (Ttrs[2],Ttrs[2]+10),
                (Ttrs[3],Ttrs[3]+20),
                (Ttrs[4],Ttrs[4]+20)]
                
                
y_plots_limits = [(0,4),(0,30),(0,60),(0,100)]

#adjust time lims and ticks 
xticks = [[xlims[1], xlims[1]+(xlims[end] - xlims[1]) / 2, xlims[end]] for xlims in x_plots_limits]

yticks = [[0,1,2,3],
    [0:10:20;],
    [0:20:50;],
    [0, 40 ,80]
    ]

for i in 1:length(plots)
    plot!(plots[i],xlims=x_plots_limits[i],
        xticks=xticks[i], yticks=yticks[i], ylims=y_plots_limits[i])
end

#adjust plots
plot!(plots[1], ylabel=L"x_i",left_margin=1Plots.mm)
plot!(plots[3], ylabel=L"x_i", left_margin=1Plots.mm)
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

#add panel level annotations 
#annotate!(plots[1], -5, y_plots_limits[1][end], text("(d)", :gray10, :right, 18))

pl4 = plot(plots...,layout=(2,2),size=(500,400), dpi = 300)

figs_dir = "../../paper_figs/"
savefig(pl4,figs_dir*"timeseries_delayed_LV_N$(N)_"*"labeled_threecolors_2x2.png")