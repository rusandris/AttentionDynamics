using DynamicalSystems
using AttentionDynamics
using Attractors
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

#find all limit cycles in sistems of size N=3,4,5
#use recurrence method from Attractors.jl


Ns = [3:5;]
diffeq = (alg=Tsit5(), dt=1e-3, adaptive=true, dtmax=1e-2, maxiters=1e15)

g_rel = 1.75 #5.0
grid_size = 1000 #phase space grid cell size
#grid_edges = [0.0,25.0]

#----------------methods parameters--------------------
T_test = 1000
Ttrs = [500,1000,1000]
cons_recurrences = 5000; 
locate_steps = 10000;
max_ensemble = 500
nr_runs = 5

#-----------critical parameters for each system---------
wc1_vals = [(N - 1) / (g_rel * (N-1)) for N in Ns]
wc2_vals = [(N - 1) / N for N in Ns]
ws = [0.6,0.6,0.7]
gs = [g_rel * (N - 1) for N in Ns]

savedir_root = "../../data/attractor_finding/"
savedir = savedir_root * "nr_attractors_cumulative_avg/"
mkpath(savedir)

#repeat attractor finding many times
attractor_finding_results = []
    
for (n,N) in enumerate(Ns)
    @show N
    g = g_rel * (N-1)
    wc1 = wc1_vals[n]
    wc2 = wc2_vals[n]
    #w = wc1 + (wc2-wc1)/2 
    w = ws[n]
    
    @show wc1,wc2,g,w
    
    #create dynamical system
    ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)
    
    #test trajectory for efficient phase space grid choice
    traj,_ = trajectory(ds,T_test;Ttr=Ttrs[n])
    #need max finding along each dim
    max_values = maximum(Matrix(traj);dims=1) #ugly max finding
    
    #construct grid in whole phase space
    grid = Tuple([range(0.0,max_values[i];length=grid_size) for i in 1:2N])
    
    #create mapper
    mapper = AttractorsViaRecurrences(ds, grid;Ttr=Ttrs[n],
    consecutive_recurrences = cons_recurrences, attractor_locate_steps = locate_steps)

    #loop over ensemble nr_runs number of times
    found_attractors = zeros(max_ensemble,nr_runs)
    for j in 1:nr_runs 
        nr_attractors = 0
        i = 0
        attr_IDs = zeros(max_ensemble)
        while (i < max_ensemble) #(nr_attractors < factorial(N-1) ) 
            attr_ID = mapper(rand(2N))
            i+=1
            attr_IDs[i] = attr_ID 
        end
        nr_attractors = length(unique(attr_IDs))
        @show nr_attractors
        found_attractors[:,j] .= [length(unique(attr_IDs[1:i])) for i in 1:max_ensemble]
    end

    #average over multiple runs
    found_attractors_avg = mean(found_attractors; dims=2)
    push!(attractor_finding_results, found_attractors_avg)
end


#save results and params 
writedlm(savedir*"found_attractors_nr_cumulative_avg_N3N4N5.txt",hcat(attractor_finding_results...))
writedlm(savedir * "systems_params.txt",hcat(ws,gs))
#writedlm(savedir*"method_params")