using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
cd(@__DIR__)

#def critical curves 
critical_curve(w, N) = (N - 1) / w
critical_curve_adim(w) = 1 / w
chaos_curve(w) = 1 / (1.5*w^2)
N = 10

#plot parameters 
ls1 = :dash
ls2 = :dot
color_default = :gray10
lw = 2
xlims = (0.0,1.0)
ylims = (0.0, 10.0)

ws = [0.1:0.01:1.0;]
gs = [1:0.1:10.0;]

#----------------plot----------------
pl = plot(guidefontsize=25, tickfontsize=18,right_margin=5Plots.mm)
plot!(pl, ws, critical_curve_adim.(ws);
    ls=ls1,lw=lw,lc=color_default, 
    label="", xlims = xlims)

vline!(pl, [(N - 1) / N],
    lc=color_default,ls=ls2,lw=lw,
    label="", xlims=xlims, ylims=ylims)

plot!(pl, ws[1:end-30], chaos_curve.(ws[1:end-30]); 
    fillrange = critical_curve_adim.(ws),
    fillstyle = :/,color=:red,
    ls=ls1, lw=lw, lc=:red,
    label="", xlims=xlims)