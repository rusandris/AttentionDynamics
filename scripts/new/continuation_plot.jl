using Revise
using BifurcationKit
using Plots,LaTeXStrings
using AttentionDynamics
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
#using DynamicalSystems
using ChaosTools
using JLD2
cd(@__DIR__)

isunity(z) = isapprox(real(z),1.0;atol=1e-1) && isapprox(imag(z),0.0;atol=1e-1)
iszero_approx(z) = isapprox(real(z),0.0;atol=1e-2) && isapprox(imag(z),0.0;atol=1e-2)

#---------------------define system-----------------------
#general so dimension can be changed

function skipsum(a,i)
	x = a[1:2:end]
    return sum(x) - a[i]
end

function F_ad(u,p)
    dim = length(u)
    N = dim/2
    du = Real[]
    for i in 1:2:(dim)
        dxi = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(N-1)*skipsum(u,i))
        dyi = u[i] - u[i+1]
        push!(du,dxi,dyi)
    end
    du
end

#-------------------helper recording and plotter functions---------------------

funcs_po = (record_from_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, p.p))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
		# plot!(br; subplot = 1, putspecialptlegend = false)
	end)

#----------define parameters------------
Ns = [3:7;] #nr of competitors dim = 2N 
cont_dir = "../../bifkit_branches/" #dir that will contain branch files 
mkpath(cont_dir)
#select g value from the phase diagram
gs_N10 = vec(readdlm("../../plotting/data/N=10.0_gvalues.txt"))
gs = gs_N10 ./(10-1)
gval_phase_diag = gs[30]
Δt = 1e-2 #lim cycle sol sampling

branches_limc = []

#loop through N values
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
    p = [g,w_max,N]

    #def problem 
    prob = BifurcationProblem(F_ad,
            u0, # initial condition x0
            p, # set of parameters
            2;     # parameter index for continuation
            )

    #------------------solve for a lim cycle---------------
    #take lim cycle closer to wc2

    #transient
    prob_de = ODEProblem(att_din!, u0, (0,5000.0), p;saveat=[5000.0])
    sol = DifferentialEquations.solve(prob_de, Rodas5())

    #higher precision, shorter traj (couple of periods)
    prob_de = ODEProblem(att_din!, sol.u[end], (0,100.), p, reltol = 1e-8, abstol = 1e-10)
    sol = DifferentialEquations.solve(prob_de, Rodas5(),saveat=Δt)
    #display(plot(sol))

    #approximate period 
    x1_ts = [sol.u[i][1] for i in 1:length(sol)]
    display(plot(sol.t,x1_ts,label="N=$N"))
    T_period = estimate_period(x1_ts,:ac,sol.t;L=length(x1_ts)÷3,ϵ=0.2)
    @show T_period

    #setup for periodic orbit (lim cycle) continuation
    #perform continuation wc_max -> w_min direction (backwards)
    probtrap, orbit_guess = generate_ci_problem(PeriodicOrbitTrapProblem(M = 250),prob, sol, T_period)
    opts_po_cont = ContinuationPar(p_min = w_min, p_max = w_max, ds = -0.02,
    newton_options = NewtonPar(;eigsolver=DefaultEig(;which=abs)), save_to_file=true,
    dsmax = 0.02,detect_bifurcation=2, nev = 2N, max_steps = 150, tol_stability = 1e-10)

    #trapezoid method continuation 
    #save max of x1
    brpo = continuation(probtrap, orbit_guess, PALC(), opts_po_cont,; 
	verbosity = 3,filename=cont_dir*"limit_cycle_branch_N$N" ,plot = false,funcs_po...)
    #finalise_solution = (z, tau, step, contResult; k...) -> begin
    #@info "Floquet exponents:"
    #(Base.display(contResult.eig[end].eigenvals) ;true)
    #end)

    push!(branches_limc,brpo)
   
end


#-------------------extract data from continuation branches-------------------------
real_parts_allN = [] #store floq multipliers real parts
imag_parts_allN = [] #store floq multipliers imag parts
sum_floq_exps_allN = [] #trace
prod_floq_exps_allN = [] #det

maxreal_floq_exps_allN = [] #maxim real part of floq exps

special_points = []
ws_step_values = []
for brpo in branches_limc
    periods = brpo.branch.period 
    @show periods[1]
    floq_exps = [brpo.eig[i].eigenvals for i in 1:length(brpo.eig)]
    floq_mult = [exp.(floqs[i].*periods[i]) for i in 1:length(floq_exps)]

    #filter out zero exponents
    floq_exps_no_zero = [filter(μ -> !(iszero_approx(μ)),floq_exps[i]) for i in 1:length(floq_exps)]
    
    #=
    #trace and det
    sum_floq_exps = sum.(floq_exps_no_zero)
    push!(sum_floq_exps_allN,sum_floq_exps)

    prod_floq_exps = prod.(floq_exps_no_zero)
    push!(prod_floq_exps_allN,prod_floq_exps)
    =#

    maxreal_floq_exps = maximum.(real.(floq_exps_no_zero))
    push!(maxreal_floq_exps_allN,maxreal_floq_exps)

    #extract w params for all continuation steps
    ws = brpo.branch.param
    push!(ws_step_values,ws)
    #should filter out unity values
    @show length(floq_mult)
    floq_mult_no_unity = [filter(λ -> !(isunity(λ)),floq_mult[i]) for i in 1:length(floq_mult)]

    real_parts = stack([real.(floq_mult[i]) for i in 1:length(floq_mult)])'
    imag_parts = stack([imag.(floq_mult[i]) for i in 1:length(floq_mult)])'

    push!(imag_parts_allN,imag_parts)
    push!(real_parts_allN,real_parts)
end

#----------------plot each branch of periodic orbits---------------------
branch_plots = []
for (i,N) in enumerate(Ns)
    #critical values
    g = gval_phase_diag * (N-1)
    wc1 = (N-1)/g
    wc2 = (N-1)/N

    #select branch object
    brpo = branches_limc[i]
    linestyle = [stb ? :solid : :dot for stb in brpo.stable]
    linewidth = [stb ? 2 : 1 for stb in brpo.stable]
    pl = plot(brpo.param, brpo.branch.max, lc=:gray10,legend=false,
        title=L"N=%$N"*"\n"*L"g/(N-1)=%$(round(gval_phase_diag;digits=2))",
        ls=linestyle,lw=linewidth,xlabel=L"w",xlims=(0.4,0.9))

    #plot some (not all) bif points
    for bp in brpo.specialpoint
        if bp.type == :ns  
            scatter!([bp.param], [brpo.branch.max[bp.idx]],ms=3,mc=:orange,markerstrokewidth=0.0)
        elseif bp.type == :nd
            scatter!([bp.param], [brpo.branch.max[bp.idx]],ms=3,mc=:blue,markerstrokewidth=0.0)
        end
    end

    #plot!(pl,brpo,linewidthunstable=0,linewidthstable=0) #plot only bif points hack

    vline!(pl,[wc1],ls=:dash,lc=:gray10)
    vline!(pl,[wc2],ls=:dash,lc=:gray10)
    i == 1 && ylabel!(pl,L"max(x_1)")

    push!(branch_plots,pl)
end
pl_bif = plot(branch_plots...,layout=(1,length(Ns)),
    size=(1500,300),left_margin=8Plots.mm,bottom_margin=8Plots.mm,top_margin=8Plots.mm)

#------------------------plot Floquet exponents/multipliers--------------------------
#plot all non-unity floquet multis (real and imag parts)    
#plot trace and det? 
#plot real part?
#plot complex plane with unit circle?
plots_imag = []
plots_real = []
plots_trdet = []
plots_maxreal = []

for (i,N) in enumerate(Ns)

    pl_maxreal = plot(ws_step_values[i],maxreal_floq_exps_allN[i],mc=:gray10,xlabel=L"w",
    ylims=(-2,2),xlims=(0.4,0.9),st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    hline!([0.0],ls=:dash)
    i == 1 && ylabel!(pl_maxreal,L"max(Re(μ_i))")

    push!(plots_maxreal,pl_maxreal)

    #=
    pl_trdet = plot(ws_step_values[i],real.(sum_floq_exps_allN[i]),mc=:gray10,xlabel=L"w",
        xlims=(0.4,0.9),st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)

    plot!(pl_trdet,ws_step_values[i],real.(prod_floq_exps_allN[i]),ls=:dash,mc=:gray10,xlabel=L"w",
        xlims=(0.4,0.9),st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    i == 1 && ylabel!(pl_trdet,L"tr,det")

    push!(plots_trdet,pl_trdet)
    
    pl_imag = plot(ws_step_values[i],imag_parts_allN[i],mc=:gray10,xlabel=L"w",
        xlims=(0.4,0.9),st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    i == 1 && ylabel!(pl_imag,L"Im(floqm_i)")
    pl_real = plot(ws_step_values[i],real_parts_allN[i],mc=:gray10,
        xlims=(0.4,0.9),st=:scatter,ms=2,markerstrokewidth=0.0,legend=false)
    push!(plots_imag,pl_imag)
    push!(plots_real,pl_real)
    =#
end
#plot(plots_real...,layout=(1,length(Ns)),size=(1500,300),ylims=[-1,1])
#pl_floq = plot(plots_imag...,layout=(1,length(Ns)),
#    size=(1200,300),ylims=[-1,1],left_margin=8Plots.mm,bottom_margin=8Plots.mm)

pl_floq = plot(plots_maxreal...,layout=(1,length(Ns)),
    size=(1200,300),left_margin=8Plots.mm,bottom_margin=8Plots.mm)
pl = plot(pl_bif,pl_floq,layout=(2,1),size=(1500,600),dpi=300)
savefig(pl,"../../continuation_test_maxreal.png")