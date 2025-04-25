using ChaosTools
using AttentionDynamics
using DynamicalSystemsBase
using OrdinaryDiffEq
using DelimitedFiles
cd(@__DIR__)
#----------define parameters------------
od_dir = "../../orbit_diagram_data/"
mkpath(od_dir)
Ns = [3:7;] #nr of competitors dim = 2N 
#select g value from the phase diagram
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
gval_phase_diag = gs[30]

n = 100 # nr of points for each param value
Ttr = 5000 #transient time 
plane = (3,1.0) #poincare plane 
comp = 1 #save components x1 of the points on the poincare surface
p_index = 2
#Δw = 0.0005 param step size
nw = 200 #param resolution in w
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
    #w = wc1 + (wc2-wc1)/2
    #@show wc1,wc2,w,g
    u0 = rand(2N)

    #params for continuation
    w_min = 1.01 * wc1
    w_max = 0.9 * wc2
    w_values = range(w_max,w_min;length=nw) 
    @show length(w_values)

    p = [g,w_max,N]


    ds = create_ds(u0,p;diffeq=diffeq) 
    pmap = PoincareMap(ds, plane)

    od = orbitdiagram(pmap, comp, p_index, w_values; n=n, Ttr=Ttr, show_progress=true)
    writedlm(od_dir*"od_N$N"*".txt",vcat(w_values',stack(od)))
    push!(orbit_diag_data_allN,od)
end