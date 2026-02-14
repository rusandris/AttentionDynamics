using DelimitedFiles
cd(@__DIR__)
include("../../scripts/new/chart_functions.jl")

#read precomputed trajectories 
#add lifetimes and save for each param 

#---------------params-----------------------


#traj specific parameters (not from filename, hardcoded)
N = 10
Δt_traj = 0.1 #sampling of trajectories
T = 100000
t = [0.0:Δt_traj:T;]
data_shape = (Int(T / Δt_traj) + 1, N)
k = 4

#trajectory path
#data_dir_traj_root = "../../data/trajs_precomp_rescaled_alpha_curve/"
data_dir_traj_root = "../../data/trajs_precomp_rescaled_vertical/"

#save path 
savedir = "../../data/lifetimes_precomp_rescaled_vertical/"
mkpath(savedir)

#read params
files = readdir(data_dir_traj_root)
params = readdlm(data_dir_traj_root * files[1])
writedlm(savedir*"B_param_pairs.txt",params)#save params also here
#params = [0.45 22.5; 0.5 22.5; 0.55 22.5]
ws = params[:,1]
gs = params[:,2]
gs_rounded = round.(gs; digits=4)
ws_rounded = round.(ws; digits=4)
chart_times = [0.1, 0.5 ,1.0]

#loop over chart times
for (idx, chart_time) in enumerate(chart_times)
    path = data_dir_traj_root 

    #---------------------------read trajectories filenames------------------------------------ 
    dir_content = readdir(path)
    files = dir_content[findall(f -> occursin(".txt",f) || occursin(".dat",f),dir_content)]

    files_trajs = files[findall(f -> occursin(".dat", f), files)] #select binary only
    nr_ics = Int(length(files_trajs)/length(gs))

    files_trajs_ics = [files_trajs[findall(f -> occursin("g$(gs[i])_w$(ws[i])", f), files_trajs)] for i in 1:length(gs)]

    #-----------------------------compute and gather lifetimes------------------------------

    #loop through parameters (corresponding filesets)
    for (i,fileset) in enumerate(files_trajs_ics)

        lifetimes_all = Int64[]

        #read each initial condition file 
        #compute lifetimes for each
        for f in fileset
            @show f
            fname, = splitext(f)
            
            #read trajectory
            traj = Matrix(reshape(reinterpret(Float64, read(path * f)), data_shape))

            #save traj pl
            mkpath(path*"figs/")
            savefig(plot(traj[1:500,:],legend=false,title=fname,size=(1200,300)),path*"figs/"*fname*".png")
            
            #integrate activity segment to segment
            int_act, time_segments = integrate_activity(traj, t, Int(chart_time / Δt_traj))
            
            lts = get_lifetimes_integ(int_act, k)
            @show length(lts)
            push!(lifetimes_all,lts...)
        end

        #save lifetimes (integers) into binary format
        base_name, = splitext(fileset[1])
        lifetimes_fname = "lifetimes_k$(k)" * base_name[6:end-3] * "charting_$chart_time" *".dat"
        write(savedir * lifetimes_fname, lifetimes_all)
    end

end