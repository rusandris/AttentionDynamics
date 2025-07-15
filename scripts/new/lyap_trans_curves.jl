using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("bif_floq_od_method_parameters.jl")


diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
#N = 10
#g = 45.0
#w = 0.66

#N = 3
#g = 10.0

N = 5
g = 20.0

#N=10 first NS bif at w = 0.329206866618409 

#w_NS = 0.329206866618409 
w_NS = 0.2599084686704565
wc1 = (N-1)/g
w_min = w_min_ratio_od * wc1
w_max = 0.4
nw = 20
ws = range(w_max,w_min;length=nw) 
u0 = rand(2N) .* 10
#u0 = generate_initial_condition(N)

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)
ds_reinit = create_ds(u0,[g,ws[1],N];diffeq=diffeq)

T = 500
Ttrs = [0,100,1000,5000,10000]

lyap_allN = [] #save lyaps here
pl = plot()
pl_reinit = plot()

for Ttr in Ttrs
    #lyapunov exps of the whole continuous system 
    lyaps = zeros(length(ws))
    lyaps_reinit = zeros(length(ws))

    #create pds 
    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9 .* rand(2N)])
    pds_reinit = ParallelDynamicalSystem(ds_reinit,[u0,u0 .+ 1e-9 .* rand(2N)])

    for (i,w) in enumerate(ws)
        @show w
        set_parameter!(pds,2,w)
        set_parameter!(pds_reinit,2,w)

        reinit!(pds_reinit)

        lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
        lyaps_reinit[i] = lyapunov(pds_reinit,T;Ttr=Ttr)
        #@show lyaps[i]
    end  
    plot!(pl,ws,lyaps,xlabel=L"w",ylabel=L"\lambda",label="Ttr=$Ttr",title="no reinit")
    plot!(pl_reinit,ws,lyaps_reinit,xlabel=L"w",ylabel=L"\lambda",label="Ttr=$Ttr",title="reinit")

end
vline!(pl,[w_NS],ls=:dash,lc=:orange,lw=2,label="NS")
vline!(pl_reinit,[w_NS],ls=:dash,lc=:orange,lw=2,label="NS")
pl_lyap = plot(pl,pl_reinit,size=(1000,400))

#some example time series, phase portraits

#ws_examples = [0.4, 0.35, 0.3, 0.25]
ws_examples = [0.22,0.25, 0.3, 0.35, 0.4]

plots_ts = []
plots_pp = []
for w in ws_examples
    set_parameter!(ds,2,w)
    traj,t = trajectory(ds,T,u0;Ttr=Ttrs[end],Δt=1e-2)

    pl_ts = plot(t[1:5000],traj[1:5000,1],lc=:gray10,title = L"w=%$w")
    pl_pp = plot(traj[:,1],traj[:,2],lc=:gray10,title = L"w=%$w")

    push!(plots_ts,pl_ts)
    push!(plots_pp,pl_pp)

end

plot_ts = plot(plots_ts...,layout=(1,length(ws_examples)),size=(1200,400))
plot_pp = plot(plots_pp...,layout=(1,length(ws_examples)),size=(1200,400))
pl_examples = plot(plots_ts...,plots_pp...,layout=(2,length(ws_examples)),size=(1200,800))


savefig(pl_lyap,"lyap_curves_N$(N)_g$g"*"_comparison.png")
savefig(pl_examples,"dynamics_examples_N$(N)_g$g"*"_comparison.png")

#=
#threaded version 

Threads.@threads for i in 1:length(ws)
    id = Threads.threadid()
    @show id
    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9])
    @show ws[i]
    set_parameter!(pds,2,ws[i])
    lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
end  
=#

#=
for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = g_rel * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #generate standard init cond
    #u0 = rand(2N)
    u0 = generate_initial_condition(N)
    @show u0

    #w_min = 1.01 * wc1
    #w_max = 0.9 * wc2
    #params for lyaps 
    w_min = w_min_ratio_od * wc1
    #w_max = w_max_ratio * wc2

    w_max = 0.4
    w_values = range(w_max,w_min;length=nw) 
    @show length(w_values)

    p = [g,w_max,N]
    #create system
    ds = create_ds(u0,p;diffeq=diffeq)

    #lyapunov exps of the whole continuous system 
    lyaps = zeros(length(w_values))

    #create pds 
    pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9])

    for (i,w) in enumerate(w_values)
        @show w
        set_parameter!(pds,2,w)
        lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
    end  

    writedlm(od_dir*"lyaps_N$N"*".txt",hcat(w_values,lyaps))
    push!(lyap_allN,lyaps)

end

=#