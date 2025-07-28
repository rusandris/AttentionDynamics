using Plots,LaTeXStrings
using DelimitedFiles
cd(@__DIR__)

data_dir = "../../data/data_highres/"

#--------------------------------------------N=5 system----------------------------------------------

#read data in
data_files_N5 = ["lyapunovs_g_20.0_N_5_Ttr_500000.txt","lyapunovs_g_20.0_N_5_Ttr_1000000.txt"]
Ttrs_strings = [L"5 \times 10^5",L"10^6"]
colors = [:gray10,:gray50]

#plot no-reinit
pl_N5_noreinit = plot(xlabel=L"w",ylabel=L"\lambda");
for (i,file) in enumerate(data_files_N5)
    data = readdlm(data_dir*file)
    params = data[:,1]
    lyaps = data[:,2]
    
    plot!(params,lyaps,lc=colors[i],label=Ttrs_strings[i])
end

#plot reinit
lyaps_N5_reinit = readdlm(data_dir*"lyapunovs_reinit_g_20.0_N_5.txt")
params = lyaps_N5_reinit[:,1]
lyaps = lyaps_N5_reinit[:,2:3] 
pl_N5_reinit = plot(params,lyaps_N5_reinit[:,2:3],
    label = [L"5 \times 10^4" L"10^5" L"2 \times 10^5"],
    lc = [:gray50 :gray30 :gray10],xlabel=L"w",ylabel=L"\lambda")

#--------------------------------------------N=10 system----------------------------------------------
data_files_N10 = ["lyapunovs_g_45.0_N_10_Ttr_50000.txt",
    "lyapunovs_g_45.0_N_10_Ttr_100000.txt",
    "lyapunovs_g_45.0_N_10_Ttr_200000.txt"]

Ttrs_strings = [L"5 \times 10^4",L"10^5",L"2 \times 10^5"]
colors = [:gray10,:gray30,:gray50]

#plot no-reinit
pl_N10_noreinit = plot(xlabel=L"w",ylabel=L"\lambda");
for (i,file) in enumerate(data_files_N10)
    data = readdlm(data_dir*file)
    params = data[:,1]
    lyaps = data[:,2]
    
    plot!(params,lyaps,lc=colors[i],label=Ttrs_strings[i])
end

#plot reinit
lyaps_N10_reinit = readdlm(data_dir*"lyapunovs_reinit_g_45.0_N_10.txt")
params = lyaps_N10_reinit[:,1]
lyaps = lyaps_N10_reinit[:,2:4] 
pl_N10_reinit = plot(params,lyaps_N10_reinit[:,2:4],
    label = [L"5 \times 10^4" L"10^5" L"2 \times 10^5"],
    lc = [:gray50 :gray30 :gray10],xlabel=L"w",ylabel=L"\lambda")

pl = plot(pl_N5_noreinit,pl_N5_reinit,pl_N10_noreinit,pl_N10_reinit;
    layout=(2,2),size=(1000,600),guidefontsize=20,
    legendfontsize=15)

figs_dir = "../../paper_figs/"
savefig(pl,figs_dir*"lyap_curves_N10_N5.png")