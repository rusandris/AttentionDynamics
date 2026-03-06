using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
cd(@__DIR__)

#match original and rescaled params 
#small and big systems 

K = 1.0;rp = 10;rc = 10;α = 0.005;c0 = 1.4;N_orig=300;
#w_B2 = 0.21 #predefine example parameter in the chaotic region

#------------------coupling strength-------------------
cs = [1.4, 2.4, 5.4] #3.0]
wg_pairs = rescaled_params.(rp,rc,K,cs,α,N_orig)
ws_horizontal = [wg_pairs[i][1] for i in 1:length(wg_pairs)]
gs_horizontal = [wg_pairs[i][2] for i in 1:length(wg_pairs)] ./ (N_orig-1)

savedir_root = "../../data/parameters/"
savedir = savedir_root * "c_article/"
mkpath(savedir)

writedlm(savedir * "params_c.txt", cs)
writedlm(savedir*"params_rp_rc_K_α_N.txt", [rp, rc, K, α, N])
writedlm(savedir*"params_w_g.txt", hcat(ws_horizontal,gs_horizontal))


#------------------small system alpha curve-------------------
#alpha curve in terms of g(w)
#same as when r changes (rp=rc=r)

alpha_curve(w, K, c) = K * c * (1 / w - 1) #g/(N-1) = 1/(ϵ*(N-1)), rp=rc

g_example = wg_pairs[1][2] ./(N-1)
_, w_B2, w_B3 = ws_horizontal
w_B1 = 0.2

K = 1.0;
#c is set so the curve goes through B2 
g_B2 = g_example
c = g_B2 / K / (1 / w_B2 - 1)
@show g_B2, alpha_curve(w_B2, K, c) #check if alpha curve passes through B2


w_B1_B3 = [w_B1:0.001:w_B3;]
gs_alpha = alpha_curve.(w_B1_B3, K, c) #alpha curve
g_B1 = gs_alpha[1]
g_B3 = gs_alpha[end]

#solve backwards for r 
N_small = 10 
α0 = 0.005
r_B1 = g_B1 * α0 * (N_orig-1)
r_B2 = g_B2 * α0 * (N_orig - 1)
rs = [g * α0 * (N_orig - 1) for g in [g_B1,g_B2,g_B3] ]

#save 
savedir_root = "../../data/parameters/"
savedir = savedir_root * "alpha_rp=rc_article/"
mkpath(savedir)

writedlm(savedir * "params_r.txt", rs)
writedlm(savedir * "params_K_c_alpha_N.txt", [K, c, α0, N])
writedlm(savedir * "params_w_g.txt", hcat([w_B1, w_B2, w_B3], [g_B1, g_B2, g_B3]))


#------------------growth rate-------------------
#rp=rc=r

K = 1.0;rp = 10;rc = 10;α = 0.005;c0 = 2.4;N=300;

rs = [9.0:1.0:12.0;]
wg_pairs = rescaled_params.(rs,rs,K,c0,α,N)
ws = [wg_pairs[i][1] for i in 1:length(wg_pairs)]
gs = [wg_pairs[i][2] for i in 1:length(wg_pairs)] ./ (N-1)


savedir_root = "../../data/parameters/"
savedir = savedir_root * "r/"
mkpath(savedir)

writedlm(savedir * "params_r.txt", rs)
writedlm(savedir*"params_K_α_c0_N.txt", [K, α, c0, N])
writedlm(savedir*"params_w_g.txt", wg_pairs)



