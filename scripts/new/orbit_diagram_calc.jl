using ChaosTools
using AttentionDynamics
using DynamicalSystemsBase
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)

#include standard parameters and init conditions
#run necessary mkpaths
include("bif_floq_od_method_parameters.jl")


#----------define local parameters------------

#dir for od and lyapunovs
od_dir = "../../orbit_diagram_data_gN5/"
mkpath(od_dir)

n = 100 # nr of points for each param value
plane = (3,1.0) #poincare plane 
comp = [1,2] #save components x1 of the points on the poincare surface
p_index = 2
#Δw = 0.0005 param step size
nw = 50 #param resolution in w
diffeq = (alg = Rodas5(),dt = 1e-3,adaptive=true,maxiters=1e15)
orbit_diag_data_allN = [] 

#-------------------orbit diagrams for each N----------------------
for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #generate standard init cond
    #u0 = rand(2N)
    u0 = generate_initial_condition(N)
    @show u0

    #params for continuation
    w_min = w_min_ratio_od * wc1
    #w_max = w_max_ratio * wc2
    w_max = 0.4
    w_values = range(w_max,w_min;length=nw) 
    @show length(w_values)

    p = [g,w_max,N]
    @show wc1,wc2,p
    #create system
    ds = create_ds(u0,p;diffeq=diffeq)
    #--------------orbit diag------------
    pmap = PoincareMap(ds, plane)

    #include u0=u0 to reset every time to u0
    #init cond is the last state from the last param
    od = orbitdiagram(pmap, comp, p_index, w_values;n=n, Ttr=Ttr, show_progress=true)  
    writedlm(od_dir*"od_N$N"*".txt",vcat(w_values',stack(od)))
    push!(orbit_diag_data_allN,od)
end

lyap_allN = [] #save lyaps here

for (i,N) in enumerate(Ns)
    @show N
    #select g for each N (system size)
    #so g/(N-1) is the same for each N
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #generate standard init cond
    #u0 = rand(2N)
    u0 = generate_initial_condition(N)
    @show u0

    #params for continuation
    w_min = 1.01 * wc1
    #w_max = 0.9 * wc2
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

