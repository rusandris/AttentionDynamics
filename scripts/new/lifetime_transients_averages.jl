using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using StatsBase
using DelimitedFiles
using Plots,LaTeXStrings
cd(@__DIR__)
include("../../plotting/defaults.jl")
include("../../plotting/plot_functions/lyapunov.jl")
include("bif_floq_od_method_parameters.jl")
include("transient_functions.jl")

diffeq = (alg = Tsit5(),dt = 1e-3,dtmax=1e-3,adaptive=true,maxiters=1e15)
N = 10
g = 45.0
#w = 0.66

#N = 3
#g = 10.0

#N = 5
#g = 20.0
ws = [0.4:-0.01:0.21;]
u0 = rand(2N)

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)

Tmax = 10000
ensemble_size = 100
lifetimes_avg = zeros(length(ws))
lyaps_avg = zeros(length(ws))

lt = simple_transient_time(ds,u0;Tmax=Tmax, threshold=1.5,
    period_tol=1e-4,max_tol=1e-2,return_ts = false)

traj,t = trajectory(ds,300,u0;Ttr=lt-100,Δt = 1e-3)
lyapunov(ds,1400;Ttr=100,u0=current_state(ds))
#traj,t = trajectory(ds,200,current_state(ds) .+ 1e-4 .* rand(2N);Ttr=0,Δt = 1e-3)

plot(t,traj[:,1])


for (j,w) in enumerate(ws)
    set_parameter!(ds,2,w)
    @show w
    lifetimes = Float64[]
    lyaps = Float64[]

    for i in 1:ensemble_size
        @show i

        #calc lifetime
        lt = simple_transient_time(ds,rand(2N);Tmax=Tmax, threshold=1.5,
            period_tol=1e-4,max_tol=Inf,return_ts = false)
        push!(lifetimes,lt)

        #calc lyap 
        #l = lyapunov(ds,lt;u0=initial_state(ds),Ttr=0)
        #push!(lyaps,l)
    end
    lifetimes_noinf = filter(x->!(isinf(x)),lifetimes) #filter Infs
    @show length(lifetimes_noinf)
    lifetimes_avg[j] = mean(lifetimes_noinf)
    #lyaps_avg[j] = mean(lyaps)
end

plot(reverse(ws),reverse(lyaps_avg),lc=:gray10)
plot(reverse(ws),reverse(lifetimes_avg),lc=:gray10,xlabel=L"w",ylabel=L"\langle \tau \rangle")
vline!([0.2],lc=:gray10,ls=:dash)

data_dir = "../../transient_lifetimes/"
mkpath(data_dir)
writedlm(data_dir*"transient_lifetime_avg_ens$ensemble_size"*"N$N"*"g_$g"*".txt",hcat(ws,lifetimes_avg))