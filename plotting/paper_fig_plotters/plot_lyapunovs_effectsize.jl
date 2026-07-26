using DelimitedFiles
using Plots, LaTeXStrings
using StatsBase
using AttentionDynamics
using LinearAlgebra
using DataFrames, CSV
cd(@__DIR__)

#path to lifetimes
#data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_Lx_small_sys/"
#data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_N50_N75_sys/"
data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_logsys_Ttr_5x1e4_0.2_0.95_w351/"
#data_dir_root = "../../data/lyapunovs_effect_size/lyapunovs_effect_size_L_Ttr_5x1e4_0.2_0.95_ws351/"
data_dirs = readdir(data_dir_root)[[1,3:11...,2]] .* "/" 

# L system constants
K = 1.0;r=10;α0 = 0.005;
p_L = [r, r, K, 2.4, α0, 300]

plots = []
#colors = distinguishable_colors(length(data_dirs)*2)[4:end]#[5:8]
#colors = [:blue :orange :green :purple]
colors = palette(:RdBu,12)[1:6]


g_rel_levels = [1.5, 2.0, 3.0, 4.0, 4.46, 5.0, 6.0, 6.69, 8.0, 9.0, 10.0] #timescale param levels on phase diag
first_zero_ws_L = Vector{Float64}[]
first_zero_ws_x = Vector{Float64}[]
lyapsL_file_idxs = 2:5
lyapsx_file_idxs = 6:9
ylabel_x = L"\Lambda_x" 
ylabel_L = L"\Lambda_{\log} / \alpha" 
ylims = (-0.1,3.0)
idx_offsets = [10,200,100,100,50,10,10,10,10,10,10]
xlimits = (0.2,1.0)

for (j,data_dir) in enumerate(data_dirs)

    #initialize fig
    pl = plot(ylabel=ylabel_L, xlabel=L"w",guidefontsize=20,tickfontsize=15,
        ylims=ylims,xlims=xlimits,
        size=(1200, 400),dpi=300)

    
    #read data files 
    files = readdir(data_dir_root * data_dir)
    files_lyaps_L = files[lyapsL_file_idxs]
    #files_lyaps_x = files[lyapsx_file_idxs]

    #@assert length(files_lyaps_L) == length(files_lyaps_x)
    nr_files = length(files_lyaps_L)

    #read params
    αs = readdlm(data_dir_root * data_dir * "gs_alphas.txt")[:,2]
    Ns = Int.(readdlm(data_dir_root * data_dir * "system_sizes.txt"))
    g_rel = g_rel_levels[j]
    first_zero_ws_level_L = Float64[] #specific to this param level
    #first_zero_ws_level_x = Float64[] #specific to this param level


    for i in 1:nr_files
        file_L = files_lyaps_L[i]
        #file_x = files_lyaps_x[i]

        @show file_L#file_x
        N = Ns[i]
        α = αs[i]

        #read files Lsystem
        lyaps_data_L = readdlm(data_dir_root * data_dir * file_L)
        lyaps_data_sorted_L = lyaps_data_L[sortperm(lyaps_data_L[:,2]),:]
        ws = lyaps_data_sorted_L[:,2]
        cs = lyaps_data_sorted_L[:,3]
        Λ_Ls = lyaps_data_sorted_L[:,end] 

        #replace blanks with NaNs (if there are)
        replace!(Λ_Ls, "" => NaN)
        Λ_Ls = Λ_Ls ./ α

        #=
        #select custom colors for some subplots
        if j == 3 
            color = :purple 
        elseif j == 8 
            color = :lime 
        else 
            color = colors[i+1]
        end
        =#
        
        color = colors[i+1]

        #plot 
        plot!(ws,Λ_Ls,lw=2,lc=color,
            title=L"y = %$(g_rel)",label="") #y = 1/(T_x(N-1))
        #plot!(ws,Λ_xs,markerstrokewidth=0.01,ms=2,st=:scatter,label="",mc=colors[i])

        #extract first zero points
        #skip zeros for low w values
        idx_offset = idx_offsets[j]
        #Lsystem
        idx_first_zero_L = findfirst(λ -> λ < 0.01, Λ_Ls[idx_offset:end]) 
        !(isnothing(idx_first_zero_L)) && push!(first_zero_ws_level_L,ws[idx_first_zero_L + idx_offset-1])

        #xsystem
        #idx_first_zero_x = findfirst(λ -> λ < 0.01, Λ_xs[idx_offset:end]) 
        #!(isnothing(idx_first_zero_x)) && push!(first_zero_ws_level_x,ws[idx_first_zero_x + idx_offset-1])

    end

    #add common legend finally
    if j == 1
        plot!(pl,legendfontsize = 15,legend=:topleft)
    end 
    push!(plots,pl)
    push!(first_zero_ws_L,first_zero_ws_level_L)
    #push!(first_zero_ws_x,first_zero_ws_level_x)

end

#add zero ws points to plots 
for (i,zpoints) in enumerate(first_zero_ws_L)
    vline!(plots[i],zpoints,ls=:dash,lc=:gray10,label="")
end

#---------------build dataframe from chaos end points and save-------------------------

chaos_ends_N20 = [first_zero_ws_L[i][1] for i in 1:length(first_zero_ws_L)] 
chaos_ends_N30 = [first_zero_ws_L[i][2] for i in 1:length(first_zero_ws_L)]
chaos_ends_N50 = [first_zero_ws_L[i][3] for i in 1:length(first_zero_ws_L)]
chaos_ends_N75 = [first_zero_ws_L[i][4] for i in 1:length(first_zero_ws_L)]


#save first zero ws values 
df = DataFrame(eps_rel = g_rel_levels,N20=chaos_ends_N20, N30=chaos_ends_N30, 
    N50=chaos_ends_N50, N75=chaos_ends_N75)

CSV.write("../../data/chaos_endpoints_lyap_curves_multiple_levels_T1e3_Ttr5x1e4_L_logsys.csv", df)


#plot all for checking purposes
pl_all = plot(plots...,layout=(2,6),size=(1800,500),dpi=300,
    left_margin=12Plots.mm,bottom_margin=12Plots.mm, right_margin=8Plots.mm)
savefig(pl_all,"../../paper_figs/preliminary_small_systems_lyap_curves_multiple_levels11_T1e3_Ttr5x1e4_L_logsys_recolor.png")

#plot only two for supplimentary fig
plots_B2 = plots[8] #top fig, higher level 
plots_B4 = plots[3] #top fig, higher level 

#edit labels, add markers
plot!(plots_B2,ylabel=L"\Lambda",ylims=(ylims[1],4.0))
plot!(plots_B4,ylabel=L"\Lambda",ylims=(ylims[1],ylims[2]+0.25))

#add markers 
#for N=20,N=30,N=50,N=75 (skip N = 10)
end_points_B2 = collect(df[8,2:end])
end_points_B4 = collect(df[3,2:end])
markersize = 8
for i in 1:length(end_points_B2)
    plot!(plots_B2,[end_points_B2[i]],[3.7],st=:scatter,ms=markersize,markershape=:circle,mc=colors[i+1],legend=:topright,label=L"N=%$(Ns[i])")
    plot!(plots_B4,[end_points_B4[i]],[3.0],st=:scatter,ms=markersize,markershape=:circle,mc=colors[i+1],legend=:topright,label="")
end

pl_lyap_curves_chaos_ends = plot(plots_B2,plots_B4,layout=(2,1),
    size=(400,600),dpi=300,right_margin=5Plots.mm)



#=
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
=#