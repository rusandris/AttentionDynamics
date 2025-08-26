using Revise
using Plots,LaTeXStrings
using AttentionDynamics
using DynamicalSystems
using LinearAlgebra

u0 = rand(10)
r=2.0;w1=0.7;w2=1.5
p = [r,w1,w2]
ds_symLV = create_ds_symLV(u0, p)

#-----------------------check eigenvalues------------------------

#fixed point values
x_fp(p, M) = 1 / (p[2] + p[3] * (M - 1))
x_1 = x_fp(p, 1)
x_N = x_fp(p, 10)

#set fixed point vecs
u0_1 = zeros(10)
u0_1[1] = x_1
u0_N = ones(10) .* x_N

#analytic jacobian
Ja_1 = build_jacobian_symLV(u0_1, p, 1)
Ja_N = build_jacobian_symLV(u0_N, p, 10)

#analytic jacobian eigvals
eigen_vals_analyticJ1 = eigvals(Ja_1)
eigen_vals_analyticJN = eigvals(Ja_N)

#check numerical Jacobian
Jf! = jacobian(ds_symLV)
J0 = zeros(10, 10)
Jf!(J0, u0_N, p)
eigen_vals_numJN = eigvals(J0)

Jf!(J0, u0_1, p)
eigen_vals_numJ1 = eigvals(J0)

#analytic eigval spectrum 
eigen_vals_analytic1 = sort(fixed_point_eigen_spectrum(10, 1, p))
eigen_vals_analyticN = sort(fixed_point_eigen_spectrum(10, 10, p))


#compare analytical and numerical eigenvalues
println("Fixed Point M=1 eigenvalues match:", all(isapprox.(eigen_vals_analytic1, eigen_vals_num1; atol=1e-2)))
println("Fixed Point M=N eigenvalues match:", all(isapprox.(eigen_vals_analyticN, eigen_vals_numN; atol=1e-2)))

plot(eigen_vals_analytic1, st=:scatter, mc=:red, ms=4, ylabel="eigenvalues", label="analytic ", title=L"M=1" * L",w_1=%$w1" * L",w_2=%$w2")
plot!(eigen_vals_num1, st=:scatter, mc=:gray10, ms=1, ylabel="eigenvalues", label="numeric ")
hline!([0],ls=:dash,lc=:gray10)

plot(eigen_vals_analyticN, st=:scatter, mc=:red, ms=4, ylabel="eigenvalues", label="analytic ", title=L"M=N"*L",w_1=%$w1"*L",w_2=%$w2" )
plot!(eigen_vals_numN, st=:scatter, mc=:gray10, ms=1, ylabel="eigenvalues", label="numeric ")
hline!([0], ls=:dash, lc=:gray10)

traj1, t = trajectory(ds_symLV, 150, u0_1 .+ randn(10) .* 1e-1; Ttr=0, Δt=1e-2)
trajN, t = trajectory(ds_symLV, 150, u0_N .+ randn(10) .* 1e-3; Ttr=0, Δt=1e-2)


plot(t,Matrix(traj1[:]))
plot(t, Matrix(trajN[:]))
hline!([x],ls=:dash,lc=:gray10)