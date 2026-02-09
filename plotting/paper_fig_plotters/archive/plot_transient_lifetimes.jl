using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

#read data
data_dir = "../../transient_lifetimes_N10_g45/transient_lifetimes/"
trs_stat = readdlm(data_dir * "transient_lifetime_avg_ens10000N10g_45.0.txt")
ws = readdlm(data_dir * "w_values.txt")

pl_lin = plot(ws,trs_stat[:,1],
    lc=:gray10,markershape=:circle,mc=:gray10,
    xlabel=L"w",ylabel=L"\langle \tau \rangle",
    guidefontsize=18,legendfontsize=10,label=L"N=10,g/(N-1)=5",
    xlims=(0.35,0.4))

pl_log = plot(ws,trs_stat[:,1], yaxis=:log,xaxis=:log,
    lc=:gray10,markershape=:circle,mc=:gray10,
    xlabel=L"w",ylabel=L"\langle \tau \rangle", xlims=(0.35,0.4),
    guidefontsize=18,legendfontsize=10,label=L"N=10,g/(N-1)=5")


plot!(ws,trs_stat[:,3], yaxis=:log,xaxis=:log,
    lc=:blue,markershape=:circle,mc=:blue,
    xlabel=L"w",ylabel=L"\langle \tau \rangle", xlims=(0.35,0.4),
    guidefontsize=18,legendfontsize=10,label=L"N=10,g/(N-1)=5")

plot!(ws,trs_stat[:,4], yaxis=:log,xaxis=:log,
    lc=:orange,markershape=:circle,mc=:orange,
    xlabel=L"w",ylabel=L"\langle \tau \rangle", xlims=(0.35,0.4),
    guidefontsize=18,legendfontsize=10,label=L"N=10,g/(N-1)=5")

pl_std = plot(ws,trs_stat[:,2], 
    lc=:gray10,markershape=:circle,mc=:gray10, yaxis=:log,xaxis=:log,
    xlabel=L"w",ylabel=L"std(\tau)", xlims=(0.35,0.4),
    guidefontsize=18,legendfontsize=10,label=L"N=10,g/(N-1)=5")

