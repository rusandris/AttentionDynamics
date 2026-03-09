using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
cd(@__DIR__)

#read precomputed trajectories 
#add lifetimes and save for each param 

#---------------params-----------------------

#traj specific parameters 
N, Δt_traj, T = readdlm(data_dir_traj * "traj_spec_params.txt")
N = Int(N)
t = [0.0:Δt_traj:T;]
data_shape = (Int(T / Δt_traj) + 1, N)
T_plot = 50

#---------------------------read trajectories filenames------------------------------------ 
dir_content = readdir(data_dir_traj)
files = dir_content[findall(f -> occursin(".txt", f) || occursin(".dat", f), dir_content)]

files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only

#search for param names
#=
file_sets = []
for j in 1:length(p_values[1])
    #put strings together
    par_names = ""
    for i in 1:length(p_values)
        par_names *= p_names[i] * "$(p_values[i][j])_"
    end
    @show par_names

    #files corresponding to a certain param set
    file_set = files_trajs[findall(f -> occursin(par_names, f), files_trajs)] 
    push!(file_sets,file_set)

end
=#
#files corresponding to a certain param set
par_names = ["B1","B2","C","B3","B4"]
file_sets = [files_trajs[findall(f -> occursin(par_names[i], f), files_trajs)] for i in 1:5] 

#loop over chart times
for (idx, chart_time) in enumerate(chart_times) 

    #-----------------------------compute and gather lifetimes------------------------------

    #loop through parameters (corresponding filesets)
    for (i,fileset) in enumerate(file_sets)

        lifetimes_all = Int64[]

        #read each initial condition file 
        #compute lifetimes for each
        for f in fileset
            @show f
            fname, = splitext(f)
            
            #read trajectory
            traj = Matrix(reshape(reinterpret(Float64, read(data_dir_traj * f)), data_shape))

            #save traj pl
            mkpath(data_dir_traj*"figs/")
            n = Int(T_plot/Δt_traj)
            savefig(plot(t[1:n],traj[1:n,:],legend=false,title=fname,size=(1200,300)),data_dir_traj*"figs/"*fname*".png")
            
            #integrate activity segment to segment
            int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt_traj))
            
            lts = get_lifetimes_integ(int_act, k)
            @show length(lts)
            push!(lifetimes_all,lts...)
        end

        #save lifetimes (integers) into binary format
        base_name, = splitext(fileset[1])
        lifetimes_fname = "lifetimes_k$(k)" * base_name[6:end-3] * "charting_$chart_time" *".dat"
        write(savedir_lt * lifetimes_fname, lifetimes_all)
    end

end

#=


#read params
files = readdir(data_dir_traj)
params = readdlm(data_dir_traj * "param_pairs.txt")
#=
#params = [0.45 22.5; 0.5 22.5; 0.55 22.5]
#params = [0.62 15.75; 0.63 15.75; 0.64 15.75]
ws = params[:,1]
gs = params[:,2]
gs_rounded = round.(gs; digits=4)
ws_rounded = round.(ws; digits=4)
=#

=#