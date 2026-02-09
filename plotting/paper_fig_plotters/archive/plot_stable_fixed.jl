using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../defaults.jl")

figs_dir = "./poster_figs/"
mkpath(figs_dir)

logocolors = Colors.JULIA_LOGO_COLORS
ms=8


function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end

N = 2
g = 5
wc2 = (N-1)/N 
wc1 = (N-1)/g 
w = 0.1

p = [g,w,N]

#u0 = generate_initial_condition(N)
#u0 = [0.8,0.0,1.3,0.0] 
#u0 = rand(2N)
u0 = [0.859425391783029,
 0.4449614962407381,
 0.30714193236605547,
 0.29700889917405093]
ds = create_ds(u0,p)



#plot fixed points
#need expressions, convenient notations

#M=1: dominance
M = 1
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
x = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u_d1 = [x,x,0,0]
u_d2 = [0,0,x,x]

#M=2: coexistence
u_c = ones(4)

#plot stable coexistence

traj,t = trajectory(ds,10;Ttr=0,Δt=0.01)
timeseries = Matrix(traj[:,[1:2:2N;]])

arrow_size = 3
ylims = (-0.1,3.0)
xlims = (-0.1,3.0)
pl = plot(timeseries[:,1],timeseries[:,2],lc=logocolors.red,
    xlabel="",ylabel=L"x_2",lw=arrow_size,
    guidefontsize=guide_fontsize,tickfontsize=tick_fontsize,
    aspect_ratio=1,
    size=(500,500),
    ylims=ylims,
    xlims=xlims,
    dpi=300,
    legendfontsize=legend_fontsize,
    bottom_margin=-2Plots.mm,
    left_margin=-2Plots.mm,
    #title="Stable coexistence \n"*L"w=%$(ds.integ.p[2]),1/\epsilon=%$g",
    label=L"w=%$w")
plot!(timeseries[80:81,1],timeseries[80:81,2],arrow=true,lw=arrow_size,lc=logocolors.red,label="")
#plot!([u_d1[1]],[u_d1[3]],st=:scatter,ms=ms,mc=logocolors.blue,markerstrokecolor=logocolors.blue,markerstrokealpha=1.0)
#plot!([u_d2[1]],[u_d2[3]],st=:scatter,ms=ms,mc=logocolors.blue,markerstrokecolor=logocolors.blue,markerstrokealpha=1.0)
plot!([u_c[1]],[u_c[3]],st=:scatter,ms=ms,mc=logocolors.red,label="")


#plot stable dominance
w=0.6
set_parameter!(ds,2,w)
u1 = [1.0,0.0,2.0,0.0]
u2 = [2.0,0.0,1.0,0.0]
#u1 = rand(4)
#u2 = rand(4)
traj1,t = trajectory(ds,10,u1;Ttr=0,Δt=0.01)
traj2,t = trajectory(ds,10,u2;Ttr=0,Δt=0.01)

l = 1000
plot!(pl,traj1[end-l:end,1],traj1[end-l:end,3],lc=logocolors.blue,
    xlabel=L"x_1",ylabel=L"x_2",
    guidefontsize=guide_fontsize,tickfontsize=tick_fontsize,
    aspect_ratio=1,
    arrow=true,
    lw=arrow_size,
    size=(500,500),
    dpi=300,
    #title="Stable dominance \n"*L"w=%$(ds.integ.p[2])",
    titlevspan=1.1,
    label=L"w=%$w")

plot!(pl,traj2[end-l:end,1],traj2[end-l:end,3],
    lc=logocolors.blue,
    label="",
    arrow=true,
    lw=arrow_size)


#plot fixed points
#need expressions, convenient notations

#M=1: dominance
M = 1
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
x = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u_d1 = [x,x,0,0]
u_d2 = [0,0,x,x]

#M=2: coexistencearrow_size
u_c = ones(4)

#plot!(timeseries[50:51,1],timeseries[50:51,2],arrow=:arrow,lw=5,lc=:gray10)
plot!([u_d1[1]],[u_d1[3]],st=:scatter,ms=ms,mc=logocolors.blue,markerstrokecolor=logocolors.blue,markerstrokealpha=1.0,label="")
plot!([u_d2[1]],[u_d2[3]],st=:scatter,ms=ms,mc=logocolors.blue,markerstrokecolor=logocolors.blue,markerstrokealpha=1.0,label="")
#plot!([u_c[1]],[u_c[3]],st=:scatter,ms=ms,mc=:white,markerstrokecolor=:gray10,markerstrokealpha=1.0)

#annotations

annotate!(1.2,1.2, text(L"u^*_c", logocolors.red, :center, 18))
annotate!(0.25,2.5, text(L"u_d^*", logocolors.blue, :center, 18))
annotate!(2.5,0.25, text(L"u_d^*", logocolors.blue, :center, 18))




#plc = plot(pl,pl2,layout=(2,1),size=(500,1000),left_margin=-10Plots.mm,bottom_margin=-6Plots.mm)

savefig(pl,figs_dir*"system_stable_fixed_N$N"*".pdf")