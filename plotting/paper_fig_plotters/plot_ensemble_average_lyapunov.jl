using Plots, LaTeXStrings
using DelimitedFiles
using StatsBase
cd(@__DIR__)

#--------------------------dirs------------------------
savedir = "../../paper_figs/"
data_dir_ens = "../../data/data_ensembles_fast_tests/"
data_dir_ben = "../../data/data_benettin/"

#--------------------------read data------------------------
files_ens = readdir(data_dir_ens)
files_ens_lyap_avg = files_ens[findall(fn -> occursin("ensemble_average",fn),files_ens)]
files_ens_lyap_avg_ordered = files_ens_lyap_avg[[1,3,5,2,4]] 

files_ben = readdir(data_dir_ben)[2]

#=
#read individual lyap exponents
files_ens_lyaps_all = files_ens[findall(fn -> occursin("lyaps_all",fn),files_ens)]
std_vals = []
for file in files_ens_lyaps_all
    lyaps_all = readdlm(data_dir_ens*file)'
    stds = std(lyaps_all;dims=1)
    push!(std_vals,stds)
end
=#

#T_init_traj_values = [100000,100000]
Ttrλ_ens_values = [100,500,1000,5000]
N = 10
g = 45.0
nr_ens = 1000

pl_lyap = plot(;dpi=300,xlabel=L"w",ylabel=L"\langle \Lambda \rangle, \Lambda",
        #title = L"N=%$N"*L" ,1/{\epsilon}=%$g"* L", T_\lambda = 10" * L", n_{ens} = %$nr_ens",
        guidefontsize=20,legendfontsize=9,
        legend_columns=2,legend=:bottomleft,
        ylims=(-0.2, 1.0), tickfontsize=10,
        xticks=[0.1:0.1:0.4;],
        xlims=(0.2,0.42));

#plot ensemble avg lyap curves
color = :gray10
las = [0.3,0.5,0.8,1.0]
mshapes = [:cross,:utriangle,:rect,:x] #:diamond]
#skip first, Ttr_lambda = 0

for (i,f) in enumerate(files_ens_lyap_avg_ordered[2:end])
   #T_init_traj = T_init_traj_values[i]
   Ttrλ_ens = Ttrλ_ens_values[i]
   data = readdlm(data_dir_ens*f)
   λs = data[:,2]
   ws = data[:,1] 
   plot!(ws,λs,lc=color,mc=color,markershape=mshapes[i],
        ms=4,markerstrokewidth=0.0,la=las[i],ma=las[i],
        label=L"T_{tr} = %$Ttrλ_ens")
end

#plot benettin curve
logocolors = Colors.JULIA_LOGO_COLORS
data = readdlm(data_dir_ben * files_ben)
ws_ben = data[:,1]
λs_ben = data[:,2]
plot!(ws_ben,λs_ben,lw=2,markershape=:circle,lc=logocolors.red,mc=logocolors.red,
        ms=4,markerstrokewidth=0.0,label=L"T_{tr}^{(\Lambda)} = 10^5")


#plot bifurcation line 
w_NS = 0.329206866618409 
vline!([w_NS],ls=:dash,lc=:gray10,lw=1,label="NS")

#savefig(pl_lyap,savedir*"lyap_ens_avg_vs_benettin_fast_tests.png")