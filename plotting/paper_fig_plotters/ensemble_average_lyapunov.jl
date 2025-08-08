using Plots, LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

savedir = "../../scripts/new/run_on_pc/"
data_dir = "../../scripts/new/run_on_pc/data/data_ensembles_fast_tests/"
ens_files = readdir(data_dir)
ens_files_Tl10 = ens_files[findall(fn -> occursin("ensemble_average",fn),ens_files)]
#ens_files_Ttrl0 = ens_files[findall(fn -> occursin("Ttrlambda_0",fn),ens_files)]
#ens_files_ordered = [ens_files_Ttrl100[2],ens_files_Ttrl0[1]]
ens_files_ordered = ens_files_Tl10[[1,3,5,2,4]] 

#T_init_traj_values = [100000,100000]
Ttrλ_ens_values = [0,100,500,1000,5000]
N = 10
g = 45.0
nr_ens = 1000

pl = plot(;dpi=300,xlabel=L"w",ylabel=L"\langle \lambda \rangle",
           title = L"N=%$N"*L" ,g=%$g"* L", T_\lambda = 10" * L", T_{i} = 10^4" * L", n_{ens} = %$nr_ens",
           guidefontsize=20,legendfontsize=8,legend_columns=-1,legend=:bottomleft,ylims = (-0.4,2.0))


for (i,f) in enumerate(ens_files_ordered)
   #T_init_traj = T_init_traj_values[i]
   Ttrλ_ens = Ttrλ_ens_values[i]
   data = readdlm(data_dir*f)
   λs = data[:,2]
   ws = data[:,1] 
   plot!(ws,λs,markershape=:circle,
          ms=2,markerstrokewidth=0.0,label=L"T_{tr}^{(\lambda)} = %$Ttrλ_ens")
end
pl
