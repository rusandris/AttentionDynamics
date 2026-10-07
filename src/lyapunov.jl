export lyapcurve,phasediag_lyapcurves
#TODO:start organizing here reusable functions from unorganized scripts

#get lyap curve for a given vertical level on phase diagram
function lyapcurve(p,nr_ws;skip=:none,savedir,T_WLC=1e3,Ttr_WLC=1e4,Ttr_fp=1e2,T_fp=10,log_sys=false,diffeq=diffeq)

    r,_,K,c_ref,α,N = p
    N = Int.(N)

    w_ref,g = rescaled_params(p...)
    g_rel = g/(N-1) #same as y

    #calc critical params
    w_c1 = 1/g_rel #HOPF 
    w_c2 = (N - 1) / N #heteroclinic

    #take w values around w_ref
    w_values = range(w_c1,w_c2;length=nr_ws)

    #set param sweep intervals
    c_start = c_convert.(w_c1,r,K,α,N)
    c_stop = c_convert.(w_values[end],r,K,α,N)
    c_values = c_convert.(w_values,r,K,α,N) #same for all system sizes

    mkpath(savedir)
    
    
    Threads.@threads for j in 1:length(c_values) 
        w = w_values[j] 
        c = c_values[j] 
        
        id = Threads.threadid() 
        @show c,id

        #define ds for each thread 
        u0_thread = rand(2N) .* 1e-2 #without the scaling, IT DOESN'T WORK
        p_new = [r, r, K, c, α, N]
        if log_sys
            u0_log = [log.(u0_thread)[1:2:end]...,u0_thread[2:2:end]...] #take log of Ls and reorder
            p_new_log = Union{Int64,Float64}[p_new[1:end-1]...,Int(N)]
            ds = create_logtransform_ds(u0_log, p_new_log;diffeq=diffeq)
        else
            ds = create_original_ds(u0_thread,p_new; diffeq=diffeq) 
        end

        #select calc time and transient time based on w
        Ttr = w_c1 < w < w_c2 ? Ttr_WLC : Ttr_fp
        Ttr_scaled = Ttr / α
        T = w_c1 < w < w_c2 ? T_WLC : T_fp
        T_scaled = T / α
        
        #skip parts of the phase diagram
        if skip==:chaos
            w > w_c1 && continue
        elseif skip==:fp
            w < w_c1 && continue
        end
        
        #calc lyapunov exponent
        #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
        @show current_parameters(ds)
        lambda = lyapunov(ds,T_scaled;Ttr=Ttr_scaled)
    
        #lock file access before writing 
        file_lock = ReentrantLock()
    
        lock(file_lock) do
            #open file here
            f = open(savedir*"lyaps_r$(r)_K$(K)_alpha$(α)_N$(N)"*".txt", "a")   
            #save lyapunov exponents
            write(f, "$(g_rel)\t$(r)\t$(w)\t$(c)\t$(lambda)\n")
            close(f) #close file here
        end
    end #of threaded param loop
end

#loop over eps params (vertical levels on phase diagram)
function phasediag_lyapcurves(Ns,g_rel_values;wstart=0.15,wstop=0.4,nr_ws=201,skip=:none,T_WLC=1e3,Ttr_WLC=1e4,Ttr_fp=1e2,T_fp=10,p_L_const = [10,10,0.005],log_sys=false,diffeq=diffeq,diffeqx=diffeqx,calc_xsys=false)
    g_rel_level_shortened = round.(g_rel_values;digits=2) 
    N = Ns[1]
    K,α0 = p_L_const

    for (gi,g_rel) in enumerate(g_rel_values)
        glevel_short = g_rel_level_shortened[gi] 
        
        #convert to r values from relative timescale levels and constant alpha 
        r = g_rel * α0 *(N-1)

      
        #calc critical params
        w_c1 = 1/g_rel #HOPF 
        w_c2 = (N - 1) / N #heteroclinic
 
        #set param sweep intervals
        w_values = range(wstart,wstop;length=nr_ws) #[w_c1:0.002:w_stop;] #[0.52:0.002:0.8;] #[w_c1:0.002:0.8;]
        c_start = c_convert.(wstart,r,K,α0,N)
        c_stop = c_convert.(w_values[end],r,K,α0,N)
        #c_values = [1.5:0.01:8.5;]
        c_values = c_convert.(w_values,r,K,α0,N) #same for all system sizes

        savedir = savedir_root * "lyapunovs_phase_diag_L_logsys_$(glevel_short)_scaled_times_dtmax_$(dtmax)/"
        mkpath(savedir)
       
        
        for i in 1:length(Ns)
            N = Ns[i]
            @info "N = $(N)"
 
            #generate standard init cond
            u0 = rand(2N)
            
            #open file here
            f = open(savedir*"lyaps_N$(N)_level_$(glevel_short)"*".txt", "a") 
            
            Threads.@threads for j in 1:length(c_values) 
                w = w_values[j] 
                c = c_values[j] 
              
                id = Threads.threadid() 
        
                #define ds for each thread 
                u0_thread = rand(2N) .* 1e-2 #without the scaling, IT DOESN'T WORK
                p_new = [r, r, K, c, α0, N]
                if log_sys
                    u0_log = [log.(u0_thread)[1:2:end]...,u0_thread[2:2:end]...] #take log of Ls and reorder
                    p_new_log = Union{Int64,Float64}[p_new[1:end-1]...,Int(N)]
                    ds = create_logtransform_ds(u0_log, p_new_log;diffeq=diffeq)
                else
                    ds = create_original_ds(u0_thread,p_new; diffeq=diffeq) 
                end
                #set_parameter!(ds, 5, c) 

                #select calc time and transient time based on w
                Ttr = w_c1 < w < w_c2 ? Ttr_WLC : Ttr_fp
                Ttr_scaled = Ttr / α0
                T = w_c1 < w < w_c2 ? T_WLC : T_fp
                T_scaled = T / α0
                
                #skip parts of the phase diagram
                if skip==:chaos
                    w > w_c1 && continue
                elseif skip==:fp
                    w < w_c1 && continue
                end
        
                @info "1/eps/(N-1) = $(glevel_short),T=$(T),Ttr=$(Ttr),r = $(r),wc1=$(w_c1),wc2=$(w_c2),w = $(w),c = $(c), thread id = $(id)" 
                flush(stderr) 
                
                #calc lyapunov exponent
                #lyaps[j] = lyapunov(pds,T;Ttr=Ttr)
                lambda = lyapunov(ds,T_scaled;Ttr=Ttr_scaled)
            
                #lock file access before writing 
                file_lock = ReentrantLock()
                
                lock(file_lock) do  
                    #save lyapunov exponents
                    write(f, "$(g_rel)\t$(r)\t$(w)\t$(c)\t$(lambda)\n")
                end
    
                if calc_xsys
                    u0_Ttrx = trajectory(dsx, 0; Ttr=Ttr)[1][1] 
                    dsx = create_ds(rand(2N), [g_rel*(N-1),w, N]; diffeq=diffeqx)
                    pdsx = ParallelDynamicalSystem(deepcopy(dsx),[u0_Ttrx,u0_Ttrx .+ d0])
                    lambdax = lyapunov(pdsx,T;Ttr=Ttr_lyap)

                     #save lyapunov exponents
                    open(savedir*"lyapsx_N$N"*".txt", "a") do f
                         write(f, "$(g_rel_calc)\t$(w)\t$(lambdax)\n")
                    end
                end #of x sys calc

            end #of threaded param loop
            close(f) #close file here
        end  #of system size
    end #of eps level param
end
