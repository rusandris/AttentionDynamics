using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
cd(@__DIR__)

#path to lifetimes
data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_Lx_small_sys/"
#data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws461_T_WLC1000_T_fp10_Ttr_WLC10000_Ttr_fp10_dtmax_0.001/"
#data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws_largerN_T_WLC1000_T_fp10_Ttr_WLC1000_Ttr_fp10_dtmax_0.001/"
#data_dir_lt = data_dir_lt_root * "lyapunovs_parallel_effect_size_ws_largerN_gnormalized_phasediag_sec1_T_WLC1000_T_fp10_Ttr_WLC1000_Ttr_fp10_dtmax_0.001/"
data_dirs = readdir(data_dir_root) .* "/"

# L system constants
K = 1.0;r=10;α0 = 0.005;
p_L = [r, r, K, 2.4, α0, 300]

plots = []
colors = distinguishable_colors(length(data_dirs)*2)
g_rel_levels = [2.0,3.0,4.46,5.0,6.69] #timescale param levels on phase diag
first_zero_ws = Vector{Float64}[]

for (j,data_dir) in enumerate(data_dirs)

    #initialize fig
    pl = plot(ylabel=L"\lambda_{L}", xlabel=L"w",guidefontsize=20,tickfontsize=15,
        legend=false,
        size=(1200, 400), bottom_margin=8Plots.mm, right_margin=8Plots.mm, 
        left_margin=8Plots.mm,dpi=300)
    
    #read data files 
    files = readdir(data_dir_root * data_dir)
    files_lyaps_L = files[2:end-6]
    αs = readdlm(data_dir_root * data_dir * "gs_alphas.txt")[:,2]
    Ns = Int.(readdlm(data_dir_root * data_dir * "system_sizes.txt"))
    g_rel = g_rel_levels[j]
    first_zero_ws_level = Float64[] #specific to this param level

    for (i,file) in enumerate(files_lyaps_L)

        N = Ns[i]
        α = αs[i]
        lyaps_data = readdlm(data_dir_root * data_dir * file)
        lyaps_data_sorted = lyaps_data[sortperm(lyaps_data[:,2]),:]
        ws = lyaps_data_sorted[:,2]
        cs = lyaps_data_sorted[:,3]
        Λ_Ls = lyaps_data_sorted[:,4] 
        
        #replace blanks with NaNs (if there are)
        replace!(Λ_Ls, "" => NaN)

        Λ_Ls = Λ_Ls ./ α

        #plot 
        plot!(ws,Λ_Ls,lw=2,label=L"N = %$(Ns[i])",title=L"\frac{1}{\varepsilon (N-1)} = %$(g_rel)")

        #extract first zero points
        #skip zeros for low w values 
        idx_offset = 10
        idx_first_zero = findfirst(λ -> λ < 0.01, Λ_Ls[idx_offset:end]) 
        !(isnothing(idx_first_zero)) && push!(first_zero_ws_level,ws[idx_first_zero + idx_offset-1])

    end

    #add common legend finally
    if j == 5
        plot!(pl,legendfontsize = 15,legend=:outertopright)
    end 
    push!(plots,pl)
    push!(first_zero_ws,first_zero_ws_level)
end

#add zero ws points to plots 
for (i,zpoints) in enumerate(first_zero_ws)
    vline!(plots[i],zpoints,ls=:dash,lc=:gray10,label="")
end

#save first zero ws values 



pl_all = plot(plots...,layout=(2,3),size=(1800,600),dpi=300)
savefig(pl_all,"../../paper_figs/preliminary_small_systems_lyap_curves_multiple_levels.png")




#example timeseries 
using OrdinaryDiffEq
using DynamicalSystems
N_test = 300

#L system
diffeq_L = (alg=Tsit5(), dt=1e-2, dtmax=1e-2, adaptive=true, maxiters=1e15)
u0_L = rand(2N_test) .* 1e-2
#fixed params 
K = 1.0;r=10;
p_L = [r, r, K, 2.4, r / (2000 / 299 * (N_test - 1)), 300]
ds_L = create_original_ds(u0_L, p_L; diffeq=diffeq_L)
traj_L, t = trajectory(ds_L, 1000, u0_L; Ttr=1000, Δt=1e-3)

pl_L = plot(t[1:25000], Matrix(traj_L)[1:25000, 1:2:end], legend=false,
    xlabel=L"t", ylabel=L"L_i", title=L"N=%$(N_test),L", guidefontsize=25)

#x system
w,g = rescaled_params(p_L...)
#init cond
η_L, η_Y = rescaling_factors(p_L...)
#u0_x = traj_L[end][1:end]
#u0_x[1:2:end] .*= 1 / η_L
#u0_x[2:2:end] .*= 1 / η_Y

p_x = [g,w, N_test]
diffeq_x = (alg=Tsit5(), dt=1e-3, dtmax=1e-3, adaptive=true, maxiters=1e15)

u0_x = u0_L
ds_x = create_ds(u0_x, p_x; diffeq=diffeq_x)
traj_x, t_x = trajectory(ds_x, 1000, u0_x; Ttr=1000, Δt=1e-3)

pl_x = plot(t_x[1:5000], Matrix(traj_x)[1:5000, 1:2:end], legend=false,
    xlabel=L"t", ylabel=L"x_i", title=L"N=%$(N_test),w=0.264", guidefontsize=25)

pl = plot(pl_L, pl_x,layout = (1,2),size=(800,400),left_margin=5Plots.mm)

#check lyapunovs 
T = 1000
Ttr = 1000
λ_L = lyapunov(ds_L, T / α; Ttr=Ttr / α)
λ_x = lyapunov(ds_x, T; Ttr=Ttr)
α = p_L[5]
λ_ratio = λ_L / λ_x

#check lyapunovs for bigger system 
p_L300 = [r, r, K, 2.4, 0.005, 300]
u0_L = rand(2*300) .* 1e-2
ds_L300 = create_original_ds(u0_L, p_L300; diffeq=diffeq_L)

w, g = rescaled_params(p_L300...)
p_x300 = [g,w,300]
ds_x300 = create_ds(u0_L, p_x300; diffeq=diffeq_x)

α0 = 0.005
T = 10
Ttr = 50
@time λ_L300 = lyapunov(ds_L300, T / α0; Ttr=Ttr/α0,Δt=1)
@time λ_x300 = lyapunov(ds_x300, T; Ttr=Ttr, Δt=1)
λ_ratio = λ_L300 / λ_x300


#check if moving to the right causes problems
cs = [2.4, 3.4, 5.4]
ws = [0.2640565204592287, 0.3370019226944242, 0.44668843025507665]

#Lsystem 
for c in cs
    set_parameter!(ds_L300,4,c)
    @time λ_L300 = lyapunov(ds_L300, T / α0; Ttr=Ttr / α0, Δt=1)
    @show λ_L300
end

#xsystem 
for w in ws
    set_parameter!(ds_x300, 2, w)
    @time λ_x300 = lyapunov(ds_x300, T ; Ttr=Ttr, Δt=1)
    @show λ_x300
end

#check different transient times 
Ttrs = [1,5,50]
set_parameter!(ds_L300,4,2.4)
n_ens = 3
for Ttr in Ttrs
    lambdas = []
    for n in 1:n_ens
        @show n 
        u0 = rand(2N_test) * 1e-2
        @time λ_L300 = lyapunov(ds_L300, T / α0; u0 = u0, Ttr=Ttr / α0, Δt=1)
        @show λ_L300
        push!(lambdas,λ_L300)
    end 
    @show mean_and_std(lambdas)    
end