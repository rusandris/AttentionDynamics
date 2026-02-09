using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
cd(@__DIR__)

N=10
u0 = rand(N)
r=1.0;w1=0.75;w2=0.3
p = [r,w1,w2]
ds_symLV = create_ds_symLV(u0, p)

#-----------------------check eigenvalues------------------------

#fixed point values
x_fp(p, M) = 1 / (p[2] + p[3] * (M - 1))
x_1 = x_fp(p, 1)
x_N = x_fp(p, N)

#set fixed point vecs
u0_1 = zeros(N)
u0_1[1] = x_1
u0_N = ones(N) .* x_N

#analytic jacobian
Ja_1 = build_jacobian_symLV(u0_1, p, 1)
Ja_N = build_jacobian_symLV(u0_N, p, N)

#analytic jacobian eigvals
eigen_vals_analyticJ1 = eigvals(Ja_1)
eigen_vals_analyticJN = eigvals(Ja_N)

#check numerical Jacobian
Jf! = jacobian(ds_symLV)
J0 = zeros(N, N)
Jf!(J0, u0_N, p)
eigen_vals_numJN = eigvals(J0)

Jf!(J0, u0_1, p)
eigen_vals_numJ1 = eigvals(J0)

#analytic eigval spectrum 
eigen_vals_analytic1 = sort(fixed_point_eigen_spectrum_symLV(N, 1, p))
eigen_vals_analyticN = sort(fixed_point_eigen_spectrum_symLV(N, N, p))


#compare analytical and numerical eigenvalues
println("Fixed Point M=1 eigenvalues match:", all(isapprox.(eigen_vals_analytic1, eigen_vals_numJ1; atol=1e-2)))
println("Fixed Point M=N eigenvalues match:", all(isapprox.(eigen_vals_analyticN, eigen_vals_numJN; atol=1e-2)))

plot(eigen_vals_analytic1, st=:scatter, mc=:red, ms=4, ylabel="eigenvalues", label="analytic ", title=L"M=1" * L",w_1=%$w1" * L",w_2=%$w2")
plot!(eigen_vals_numJ1, st=:scatter, mc=:gray10, ms=1, ylabel="eigenvalues", label="numeric ")
hline!([0],ls=:dash,lc=:gray10)

plot(eigen_vals_analyticN, st=:scatter, mc=:red, ms=4, ylabel="eigenvalues", label="analytic ", title=L"M=N"*L",w_1=%$w1"*L",w_2=%$w2" )
plot!(eigen_vals_numJN, st=:scatter, mc=:gray10, ms=1, ylabel="eigenvalues", label="numeric ")
hline!([0], ls=:dash, lc=:gray10)


#----------------------------symmetric no delay LV------------------------------

#WTA

#calc short trajectory
set_parameter!(ds_symLV, 2, w2)
set_parameter!(ds_symLV, 3, w1)
u0 = rand(N)*3
x_1 = x_fp(current_parameters(ds_symLV), 1)
traj1, t = trajectory(ds_symLV, 100, u0; Ttr=0, Δt=1e-2)

#set colors (first three)
idxs = sortperm(u0)[end-2:end]

colors = collect(Colors.JULIA_LOGO_COLORS)[[1,2,4]]

pl_1 = plot(xlabel=L"t",
    guidefontsize=18, legendfontsize=9,
    tickfontsize=10, xlims=(0, 50), 
    ylims=(0.0, 6.0), dpi=300)
    #title=L"N = %$N" * L", w_1 = %$w1" * L", w_2 = %$(ws[2])")

#plot top 3 with colors 
plot!(t, Matrix(traj1[:,idxs]), lw=4, la=1, label="",lc=colors')
#plot other with gray
idxs_others = []
for i in 1:N
    i ∈ idxs && continue
    push!(idxs_others, i)
end
plot!(t, Matrix(traj1[:, idxs_others]), lw=1, la=1, label="", lc=:gray50)

hline!([x_1], ls=:dash, lc=:gray10, lw=2, label="")
#L"x^*_N = \frac{1}{w_1 + (N-1)w_2}"



#coexistence
ws = [0.3,0.8]

set_parameter!(ds_symLV, 2, w1)
set_parameter!(ds_symLV, 3, w2)
x_N = x_fp(current_parameters(ds_symLV), N)

trajN, t = trajectory(ds_symLV, 100, rand(10); Ttr=0, Δt=1e-2)

pl_N = plot(ylabel=L"x_i", xlabel=L"t",
    guidefontsize=18, legendfontsize=9, 
    tickfontsize=10,xlims=(0,50),
    ylims=(0.0, 1.0), dpi=300)
    #title=L"N = %$N"*L", w_1 = %$w1"*L", w_2 = %$(ws[1])")

#plot top 3 with colors 
plot!(t, Matrix(trajN[:, idxs]), lw=4, label="", lc=colors')

#plot other with gray
plot!(t, Matrix(trajN[:, idxs_others]), lw=1, label="", lc=:gray50)

hline!([x_N], ls=:dash, lc=:gray10,lw=2, label="")
#L"x^*_N = \frac{1}{w_1 + (N-1)w_2}")

#save figs 
pl = plot(pl_N,pl_1,size=(600,300), dpi = 300,
    left_margin=2Plots.mm,bottom_margin=5Plots.mm)

figs_dir = "../../paper_figs/"

savefig(pl_N, figs_dir * "timeseries_symmetric_LV_N$N" * "_coexistence.png")
savefig(pl_1, figs_dir * "timeseries_symmetric_LV_N$N" * "_WTA.png")

savefig(pl,figs_dir*"timeseries_symmetric_LV_N$N"*".png")


#-------------------------------------phase diagram---------------------------------
#try w1 and w2 phase diagram 
w_range = [0.0:1e-2:1.0;]

phase_diag = plot(w_range, w_range, fillrange=[0.0], fillstyle=:\,
    fillcolor=:gray50,
    xlabel=L"w_1", ylabel=L"w_2",
    xlims=(0.0, 1.0), ylims=(0.0, 1.0),
    ticks = [0.0,0.5,1.0],lw=5,lc=:gray10,
    legend=false,guidefontsize=18,tickfontsize=10,
    aspect_ratio=1.0, size=(600,600),dpi=300)


#plot parameter values 

plot!(phase_diag, [w1], [w2],st=:scatter,
    mc=:white,markershape=:circle,ms=8)
plot!(phase_diag, [w2], [w1], st=:scatter, 
    mc=:gray10, markershape=:circle, ms=8)

annotate!(w1, w2+0.1, text("A", :gray10, :center, 20))
annotate!(w2, w1+0.1, text("B", :gray10, :center, 20))

savefig(phase_diag, figs_dir * "phase_diag_symmetric_LV_N$N" * ".png")
