using DelimitedFiles
using Plots,LaTeXStrings
using AttentionDynamics
cd(@__DIR__)

#match original and rescaled params 
#small and big systems 

K = 1.0;rp = 10;rc = 10;α = 0.005;c0 = 1.4;N=300;

#------------------coupling strength-------------------
cs = [1.4, 2.2, 3.0]
wg_pairs = rescaled_params.(rp,rc,K,cs,α,N)

savedir_root = "../../data/parameters/"
savedir = savedir_root * "c/"
mkpath(savedir)

writedlm(savedir * "params_c.txt", cs)
writedlm(savedir*"params_rp_rc_K_α_N.txt", [rp, rc, K, α, N])
writedlm(savedir*"params_w_g.txt", wg_pairs)

#------------------growth rate-------------------
#rp=rc=r

K = 1.0;rp = 10;rc = 10;α = 0.005;c0 = 2.4;N=300;

rs = [9.0:1.0:12.0;]
wg_pairs = rescaled_params.(rs,rs,K,c0,α,N)

savedir_root = "../../data/parameters/"
savedir = savedir_root * "r/"
mkpath(savedir)

writedlm(savedir * "params_r.txt", rs)
writedlm(savedir*"params_K_α_c0_N.txt", [K, α, c0, N])
writedlm(savedir*"params_w_g.txt", wg_pairs)
