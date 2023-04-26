using AttentionDynamics,OrdinaryDiffEq,DelimitedFiles

#-------------------system size and parameters-----------------------
N = 4
g = 8.0
wstart = 0.68
wstop = 0.71
p = [g,wstart,N]

#----------------------select initial condition----------------------
u0s = readdlm("N4_LIMCYCLES_IC.txt")
ic_index = 5
u0 = u0s[:,ic_index]


#poincare plane x2 == 1.0 
plane = (3,1.0)

#pvalues from w_start,w_stop and Δw
Δw = 0.002
pvalues = [wstart:Δw:wstop;]

#-----------------------comp and transient time----------------------
Ttr = 1000
tfinal = 1000

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=false,maxiters=1e15)

outfilename = "od_data_N$N"*"g_$g"*"_w_$wstart"*"_$wstop"*"_tf$tfinal"*"_Ttr$Ttr"*"_ic$ic_index"*".txt"

#---------------------orbit diagram calculation------------------------
ds = create_ds(u0,p)
od = orbit_diagram(ds,plane,1,2,pvalues;direction = 1, Ttr = Ttr,tf = tfinal,diffeq=diffeq,filename=outfilename)
save_od_data(outfilename,od)
