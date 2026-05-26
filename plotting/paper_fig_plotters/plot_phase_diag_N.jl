using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
cd(@__DIR__)

#def critical curves 
critical_curve(w, N) = (N - 1) / w
critical_curve_adim(w) = 1 / w
chaos_curve(w, beta) = w^(-beta)

chaos_points1 = [0.45, 0.55, 0.6, 0.52, 0.6]
rescaled_time1 = 6.88 * ones(length(chaos_points1))
Ns1 = [30, 50, 75, 200, 300]

chaos_points2 = [0.39, 0.4, 0.45, 0.52]
rescaled_time2 = 4.0 * ones(length(chaos_points2))
Ns2 = [8, 10, 15, 20]


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
pl = plot(guidefontsize=25, tickfontsize=18,right_margin=5Plots.mm,legendfontsize=12,dpi=300)
plot!(pl, ws, critical_curve_adim.(ws);
    ls=ls1,lw=lw,lc=color_default, 
    label="", xlims = xlims)

#read chaos end points for different system sizes 
chaos_end_points = readdlm("../../data/small_systems_lyap_drops_first_zeros_N_grel_w_T1e3_Ttr1e4.txt")

#go through chaos points 
colors = [:blue :orange :green :purple ]
nr_points = size(chaos_end_points)[1]
for i in 1:4:nr_points
    points = chaos_end_points[i:i+3,2:3]
    plot!(points[:,2]',points[:,1]',st=:scatter,ms=6,label="",mc=colors)
end
pl

vline!(pl, [(N - 1) / N],
    lc=color_default,ls=ls2,lw=lw,
    label="", xlims=xlims, ylims=ylims)

Ns = [20,30,50,75]
βs = [2.0,2.4,2.8,4.0]
for (i,β) in enumerate(βs)
    plot!(pl, ws, chaos_curve.(ws, β); 
        ls=ls1, lw=lw, lc=colors[i],
        label=L"N = %$(Ns[i]),\beta=%$(β)")
end
pl
#------------------------add manually----------------------

scatter!(pl, chaos_points1, rescaled_time1, color=:red, label="")
scatter!(pl, chaos_points2, rescaled_time2, color=:orange, label="")


plot!(pl, ws, chaos_curve.(ws, 1.5); 
    fillrange = critical_curve_adim.(ws),
    fillstyle = :/,color=:red,
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=1.5")
    
plot!(pl, ws, chaos_curve.(ws, 1.75); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=1.75")
    

plot!(pl, ws, chaos_curve.(ws, 2.1); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=2.1")


plot!(pl, ws, chaos_curve.(ws, 2.4); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=2.4")


plot!(pl, ws, chaos_curve.(ws, 2.95); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=2.95")

plot!(pl, ws, chaos_curve.(ws, 3.2); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=3.2")

plot!(pl, ws, chaos_curve.(ws, 3.8); 
    ls=ls1, lw=lw, lc=:red,
    label=L"\beta=3.8")

    
plot!(pl, xlims=xlims, xlabel = L"w", ylabel = L"T_x/(N-1)")


#----------------plot----------------
betas = [1.5, 1.75, 2.1, 2.4, 2.95, 3.2, 3.8]

pl2 = plot(guidefontsize=25, tickfontsize=18,right_margin=5Plots.mm)
scatter!(pl2, cat(Ns2,Ns1; dims=1)[1:end-2], betas, label="")
plot!(pl2, xscale=:log, xlabel=L"N", ylabel=L"\beta")
