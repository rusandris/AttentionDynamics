using Plots, LaTeXStrings
#using AttentionDynamics
using DelimitedFiles
using DataFrames, CSV
using StatsBase
using CurveFit  
cd(@__DIR__)
include("defaults.jl")

#def critical curves 
critical_curve(w, N) = (N - 1) / w
critical_curve_adim(w) = 1 / w
chaos_curve(w, beta, gamma=0.0) = w^(-beta) + gamma
powerlaw_translated(p, x) = @. x^(p[1]) + p[2] 
powerlaw(p, x) = @. p[2]*x^(p[1])  


#plot parameters 
ls1 = :dash
ls2 = :dot
color_default = :gray10
colors = palette(:RdBu,12)[1:6]


opac_alphas = [1.0,1.0,1.0,1.0]
#colors = [:blue :orange :green :purple ]
markershapes = [:utriangle, :rect, :cross, :circle]
#markershapes = [:circle, :circle, :circle, :circle ]
lw = 3
ms = 4
xlims = (0.0,1.0)
ylims = (0.0, 10.5)

ws = [0.1:0.01:1.0;]
gs = [1:0.1:10.0;]

#----------------------------plot----------------------
ylabel = L"y" #L"1/(T_x(N-1))" #L"\frac{1}{\varepsilon(N-1)}"
pl_phasediag_N = plot(ylabel = ylabel,xlabel=L"w",guidefontsize=20, 
    tickfontsize=12,left_margin=5Plots.mm, right_margin=5Plots.mm, 
    top_margin=8Plots.mm,legendfontsize=10, dpi=300, legend=:topright,
    legendcolumns=2)

plot!(pl_phasediag_N, ws, critical_curve_adim.(ws);
    ls=ls1,lw=lw,lc=color_default, 
    label="", xlims = xlims,ylims=ylims)

#read chaos end points for different system sizes 
chaos_end_points = CSV.read("../../data/lyapunovs_effect_size/chaos_endpoints_lyap_curves_multiple_levels_T1e3_Ttr5x1e4_L_logsys.csv",DataFrame)
nr_points = size(chaos_end_points)[1]

#go through chaos points 
g_rel_levels = chaos_end_points[:,1]
Ns = [20,30,50,75] #unique(chaos_end_points[:,1])

#add horizontal levels 
hline!([g_rel_levels[8]],ls=:dash,lw=3,lc=:lime,label="")#label=L"6.68")
hline!([2.98],ls=:dash,lw=3,lc=:purple,label="")#label=L"2.98")


#read phase diagram of N=10 separately
data_dir_N10 = "../../data/phasediag_data_N10_extended5_80x80/"
lyap_vals_N10 = reshape(readdlm(data_dir_N10 * "lyap_values_Ttr10000.txt"), (80, 80))
params = readdlm(data_dir_N10*"phasediag_params.txt")
ws_N10 = params[:,1]
gs_N10 = params[:, 2]

#loop through every row on the phase diagram, search for end of chaos 
chaos_end_idxs_N10 = zeros(Int64,length(ws_N10))
chaos_threshold = 0.1
periodic_threshold = 1e-2
range_avg = 3
for (i,g) in enumerate(gs_N10)
    lyap_curve = lyap_vals_N10[i,:]

    #search for zero exponent after a positive one
    for j in range_avg:length(lyap_curve)-range_avg
        avg_before = mean(lyap_curve[(j-range_avg+1):j])
        avg_after = mean(lyap_curve[j:j+range_avg-1])
        if  avg_before > chaos_threshold && 0 < abs(avg_after) < periodic_threshold        
            chaos_end_idxs_N10[i] = j
        end
    end
end
#find first nonzero chaos end point 
first_nonzero_idx = findfirst(i -> !(iszero(i)),chaos_end_idxs_N10) 
chaos_end_points_N10 = ws_N10[chaos_end_idxs_N10[first_nonzero_idx:end]]

#plot N=10 separately
g_rel_levels_N10 = gs_N10[end-length(chaos_end_points_N10)+1:end] ./ (10-1)
plot!(pl_phasediag_N,chaos_end_points_N10[1:5:end],g_rel_levels_N10[1:5:end],
    st=:scatter,markershape=:circle,ms=ms,mc=colors[1],label=L"N = 10",alpha=1.0)

#add N=10 fit 
param_guesses = [-1.0, 0.5]
prob_N10 = NonlinearCurveFitProblem(powerlaw, param_guesses, chaos_end_points_N10, g_rel_levels_N10)
ps_fit_N10 = solve(prob_N10).u
β_N10, a_N10 = ps_fit_N10
println("Fit params, N = 10, beta = $(β_N10), a = $(a_N10)")
plot!(ws_N10, powerlaw(ps_fit_N10, ws_N10), lw=lw, label="", lc=colors[1], alpha=1.0)


#fit to each N 
param_guesses = [-1.0, 0.5]
params_fit = []
for (i,N) in enumerate(Ns)
    w_chaos_ends = chaos_end_points[:,i+1]
    prob = NonlinearCurveFitProblem(powerlaw, param_guesses, w_chaos_ends, g_rel_levels)
    sol = solve(prob)
    ps_fit = sol.u
    β,a = ps_fit   
    println("Fit params, N = $(N), beta = $(β), a = $(a)")
    push!(params_fit,ps_fit)
end

#plot_chaos points 
for i in 1:length(Ns)
    chaos_end_points_N = chaos_end_points[:,i+1]
    p_fit = params_fit[i]
    #strokewidth = i == 3 ? 4 : 0.0 
    plot!(chaos_end_points_N,g_rel_levels,st=:scatter,ms=ms,alpha=opac_alphas[i],
        label=L"N = %$(Ns[i])",mc=colors[i+1],markershape=:circle)
    plot!(ws,powerlaw(p_fit,ws),lw=lw,label="",lc=colors[i],alpha=opac_alphas[i])
end

#add annotations
annotate!(pl_phasediag_N, 0.2, 1.2, text("coexistence", :gray10, :center, annot_size-10))
annotate!(pl_phasediag_N, 0.32, 4.2, text("chaos", :gray10, :center, annot_size-10))
annotate!(pl_phasediag_N, 0.3, 11.0, text("chaotic WLC boundaries", :gray10, :center, annot_size-10))


savefig(pl_phasediag_N,"../../paper_figs/phase_diag_size_effect_chaos_endpoints.png")

#add heteroclinic bifurcation line 
#vline!(pl_phasediag_N, [(N - 1) / N],
#    lc=color_default,ls=ls2,lw=lw,
#    label="", xlims=xlims, ylims=ylims)

#savefig(pl_phasediag_N, "../../paper_figs/chaos_endpoints_lyap_curves_4system_sizes_multiple_levels_T1e3_Ttr5x1e4_L_logsys.png")

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
=#