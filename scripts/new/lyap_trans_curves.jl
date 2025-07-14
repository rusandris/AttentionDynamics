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

wc1 = (N-1)/g
w_min = w_min_ratio_od * wc1
w_max = 0.4
nw = 20
ws = range(w_max,w_min;length=nw) 
u0 = rand(2N) .* 10

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)

T = 2000
Ttr = 50000

lyap_allN = [] #save lyaps here


#lyapunov exps of the whole continuous system 
lyaps = zeros(length(ws))

#create pds 
pds = ParallelDynamicalSystem(ds,[u0,u0 .+ 1e-9])

for (i,w) in enumerate(ws)
    @show w
    set_parameter!(pds,2,w)
    lyaps[i] = lyapunov(pds,T;Ttr=Ttr)
end  
#plot(ws,lyaps,xlabel=L"w",ylabel=L"\lambda")

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