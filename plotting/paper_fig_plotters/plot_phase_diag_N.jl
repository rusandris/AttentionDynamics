using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
using CurveFit  
cd(@__DIR__)

#def critical curves 
critical_curve(w, N) = (N - 1) / w
critical_curve_adim(w) = 1 / w
chaos_curve(w, beta, gamma=0.0) = w^(-beta) + gamma

powerlaw(p, x) = @. x^(p[1]) + p[2] 

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
ylims = (0.0, 10.5)

ws = [0.1:0.01:1.0;]
gs = [1:0.1:10.0;]

#----------------plot----------------
pl = plot(ylabel = L"\frac{1}{\varepsilon(N-1)}",xlabel=L"w",guidefontsize=25, tickfontsize=18,
    left_margin=5Plots.mm, right_margin=5Plots.mm, legendfontsize=12, dpi=300)
plot!(pl, ws, critical_curve_adim.(ws);
    ls=ls1,lw=lw,lc=color_default, 
    label="", xlims = xlims,ylims=ylims,legend=:bottomleft)

#read chaos end points for different system sizes 
chaos_end_points = CSV.read("../../data/chaos_endpoints_lyap_curves_multiple_levels_T1e3_Ttr5x1e4_L_logsys.csv",DataFrame)
nr_points = size(chaos_end_points)[1]

#go through chaos points 
g_rel_levels = chaos_end_points[:,1]
Ns = [20,30,50,75] #unique(chaos_end_points[:,1])
#chaos_end_points_per_size = [chaos_end_points[i:4:end,end] for i in 1:4]

#fit to each N 
param_guesses = [-1.0, 0.5]
params_fit = []
for (i,N) in enumerate(Ns)
    w_chaos_ends = chaos_end_points[:,i+1]
    prob = NonlinearCurveFitProblem(powerlaw, param_guesses, w_chaos_ends, g_rel_levels)
    sol = solve(prob)
    ps_fit = sol.u
    push!(params_fit,ps_fit)
end

#plot_chaos points 
colors = [:blue :orange :green :purple ]
for i in 1:length(Ns)
    chaos_end_points_N = chaos_end_points[:,i+1]
    p_fit = params_fit[i]
    plot!(chaos_end_points_N,g_rel_levels,st=:scatter,ms=6,label="",mc=colors[i])
    plot!(ws,powerlaw(p_fit,ws),lw=1,label="",lc=colors[i])
end
pl

#hline!([1.5],ls=:dash)
vline!(pl, [(N - 1) / N],
    lc=color_default,ls=ls2,lw=lw,
    label="", xlims=xlims, ylims=ylims)

savefig(pl, "../../paper_figs/chaos_endpoints_lyap_curves_4system_sizes_multiple_levels_T1e3_Ttr5x1e4_L_logsys.png")

#=
Ns = [20,30,50,75]
βs = [1.8,2.0,2.8,3.0]
γs = [0.3,0.2,1.2,1.5]
for (i,β) in enumerate(βs)
    plot!(pl, ws, chaos_curve.(ws, β,γs[i]); 
        ls=ls1, lw=lw, lc=colors[i],
        label=L"N = %$(Ns[i]),\beta=%$(β)")
end
pl
#add manually some (test grel8)
plot!([0.32],[8.0],ms=6,mc=colors[1],st=:scatter,label="")
plot!([0.38],[8.0],ms=6,mc=colors[2],st=:scatter,label="")
plot!([0.45],[8.0],ms=6,mc=colors[3],st=:scatter,label="")
plot!([0.53],[8.0],ms=6,mc=colors[4],st=:scatter,label="")
=#

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
