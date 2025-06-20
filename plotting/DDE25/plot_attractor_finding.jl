using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../defaults.jl")
logocolors = Colors.JULIA_LOGO_COLORS
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
figs_dir = "./poster_figs/"
Ns = [3:5;]
#gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
#gs = gs_N10 ./(10-1)
#g_N = gs[30]
g_rel = 5.0

w = 0.5
grid_size = 1000 #phase space grid cell size
grid_edges = [0.0,25.0]

T_test = 1000
Ttr=5000
cons_recurrences = 5000; 
locate_steps = 10000;
max_ensemble = 100
nr_runs = 5

#=
grid = (
    range(0.0, 5.0; length = 3), # x1
    range(0.0, 5.0; length = 3), # x2
    range(0.0, 5.0; length = 3), # x3
    range(0.0, 3.0; length = 3), # y1
    range(0.0, 3.0; length = 3), # y2
    range(0.0, 3.0; length = 3), # y3
)
=#

#repeat attractor finding many times

#=
attractor_finding_results = []
    
for (n,N) in enumerate(Ns)
    @show N
    g = g_rel * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N
    w = wc1 + (wc2-wc1)/2 
    
    @show wc1,wc2,g,w
    
    ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)

    
    #test trajectory for efficient phase space grid choice
    traj,_ = trajectory(ds,T_test;Ttr=Ttr)
    #need max finding along each dim
    max_values = maximum(Matrix(traj);dims=1) #ugly max finding
    
    grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])
    
    
    mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttr,
    consecutive_recurrences = cons_recurrences, attractor_locate_steps = locate_steps)

    found_attractors = []
    for j in 1:nr_runs 
        nr_attractors = 0
        i=0
        attr_IDs = Int[]
        while (i < max_ensemble) #(nr_attractors < factorial(N-1) ) 
            attr_ID = mapper(rand(2N))
            nr_attractors = length(mapper.bsn_nfo.attractors)
            i+=1
            @show attr_ID
            push!(attr_IDs,attr_ID)
        end
        push!(found_attractors,attr_IDs)
    end
    push!(attractor_finding_results,found_attractors)
end
=#

#read attractor finding results
nr_atrs_avgs = readdlm("attractor_finding/attractor_finding_Ns345.txt")
method_params = readdlm("attractor_finding/attractor_finding_method_params.txt")
Ns = Int.(method_params[1:3])

#--------------------------------plot attractor finding----------------------------
lcolors = [logocolors.red,logocolors.green,logocolors.purple]
pl_atr = plot(;ylabel="# attractors",xlabel="# initial conditions",
    legend=:bottomright,legend_columns=-1,DEFAULTS...,legendfontsize=legend_fontsize+5)


#nr_atrs_avgs = []
#for (i,results) in enumerate(attractor_finding_results)
for (i,N) in enumerate(Ns)

    #=
    #take average of many runs
    nr_atrs_avg = zeros(max_ensemble)
    for j in 1:length(results)
        atrs = results[j]
        nr_atrs = [length(unique(atrs[1:i])) for i in 1:length(atrs)]
        nr_atrs_avg .+= nr_atrs
    end
    push!(nr_atrs_avgs,nr_atrs_avg ./ length(results))
    =#

    nr_atrs_avg = nr_atrs_avgs[:,i]
    #------------------plot----------------------------
    nr_atr_max = factorial(Ns[i]-1)
    annotate!(70, nr_atr_max*1.2, text(L"(N-1)! = %$nr_atr_max", lcolors[i], :center, 28))
    #plot!(nr_atrs_avg ./ length(results),lw=4,lc=lcolors[i],label=L"N=%$(Ns[i])")
    plot!(nr_atrs_avg,lw=6,lc=lcolors[i],label=L"N=%$(Ns[i])")
    hline!([nr_atr_max],lw=3,ls=:dash,lc=lcolors[i],label="")
end
ylims!(pl_atr,(1,30),yticks=[1,10^0.5,10],yaxis=:log,tickfontsize=tick_fontsize_large+5,guidefontsize=guide_fontsize_large+5)
plot!(size=(800,700))


savefig(pl_atr,figs_dir*"system_multistability_nr_attractors.png")
#=
data_dir = "./attractor_finding/" 
mkpath(data_dir)
writedlm(data_dir*"attractor_finding_Ns345.txt",stack(nr_atrs_avgs))
writedlm(data_dir*"attractor_finding_method_params.txt",[Ns...,g_rel,w,cons_recurrences,locate_steps,max_ensemble,nr_runs])
=#