using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra
cd(@__DIR__)

N=10
u0 = rand(N)
r=2.0;w1=0.5;w2=0.3
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


#--------------------------------plot time series---------------------------
ws = [0.3,0.8]

set_parameter!(ds_symLV, 3, ws[1])
x_N = x_fp(current_parameters(ds_symLV), N)

trajN, t = trajectory(ds_symLV, 100, rand(N); Ttr=0, Δt=1e-2)

pl_N = plot(ylabel=L"x_i", xlabel=L"t",
    guidefontsize=18, legendfontsize=18, tickfontsize=10,
    title=L"N = %$N"*L", w_1 = %$w1"*L", w_2 = %$(ws[1])")

hline!([x_N], ls=:dash, lc=:gray10, label=L"u_N")
plot!(t, Matrix(trajN[:]), label="")

set_parameter!(ds_symLV, 3, ws[2])
x_1 = x_fp(current_parameters(ds_symLV), 1)

traj1, t = trajectory(ds_symLV, 100, rand(N); Ttr=0, Δt=1e-2)


pl_1 = plot(ylabel=L"x_i", xlabel=L"t",legend=:bottomright,
    guidefontsize=18, legendfontsize=18, tickfontsize=10,
    title=L"N = %$N"*L", w_1 = %$w1"*L", w_2 = %$(ws[2])")

hline!([x_1], ls=:dash, lc=:gray10, label=L"u_1")
plot!(t, Matrix(traj1[:]), ylabel=L"x_i", xlabel=L"t",
    guidefontsize=18, legendfontsize=18, tickfontsize=10, label="")


pl = plot(pl_N,pl_1,size=(800,300), dpi = 300,
    left_margin=5Plots.mm,bottom_margin=5Plots.mm)

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"timeseries_symmetric_LV_N$N"*".png")