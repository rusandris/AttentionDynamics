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

#N=10 first NS bif at w = 0.329206866618409 
ws = [0.4:-0.001:0.362;]
u0 = rand(2N)

ds = create_ds(u0,[g,ws[1],N];diffeq=diffeq)

Tmax = 20000
ensemble_size = 10000
threshold = 1.5

fig_dir = "../../lifetime_distribs/"
mkpath(fig_dir)

avgs = zeros(length(ws))
stds = similar(avgs)
nr_infs = similar(avgs)
mins = similar(avgs)
maxs = similar(avgs)

data_dir = "../../transient_lifetimes/"
mkpath(data_dir)

file_name = data_dir*"transient_lifetime_avg_ens$ensemble_size"*"N$N"*"g_$g"*".txt"

for (j,w) in enumerate(ws)
    set_parameter!(ds,2,w)
    @show w
    init_states = [rand(20) .* 10 for i in 1:ensemble_size]

    lts = transient_lifetime_statistics(ds,init_states;
        Tmax=Tmax, period_tol=1e-3, max_tol=1e-3)

    #plot distribution
    pl = histogram(lts;bins=100,xlims=(0,6000),ylabel="",xlabel=L"\tau")
    savefig(pl,fig_dir*"lifetime_distrib_w$w"*".png")


    #save lifetime statistics
    avgs[j] = mean(lts) 
    stds[j] = std(lts)
    nr_infs[j] = count(x -> isinf(x),lts)
    min,max = extrema(lts)
    mins[j] = min
    maxs[j] = max

    open(file_name, "a") do f
        write(f, "$(avgs[j]) \t $(stds[j]) \t $(mins[j]) \t $(maxs[j]) \t $(nr_infs[j]) \n")
    end
end

#plot(reverse(ws),reverse(lyaps_avg),lc=:gray10)
#plot(reverse(ws),reverse(lifetimes_avg),lc=:gray10,xlabel=L"w",ylabel=L"\langle \tau \rangle")
#vline!([0.2],lc=:gray10,ls=:dash)


