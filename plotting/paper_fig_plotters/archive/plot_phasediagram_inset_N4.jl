using Plots,AttentionDynamics,DelimitedFiles,LaTeXStrings

include("defaults.jl")
include("plot_functions/lyapunov.jl")


output_filename = "figs/phase_diagram_N4_T=1000_Ttr=10000_res_40x40.pdf"
Λ = readdlm("data/N=4.0_T=1000_Ttr=10000lyaps_res_40x40.txt")
wvalues = vec(readdlm("data/N=4.0_wvalues_T=1000_Ttr=10000_res_40x40.txt"))
gvalues = vec(readdlm("data/N=4.0_gvalues_T=1000_Ttr=10000_res_40x40.txt"))
ns_data = readdlm("data/ns_points.txt") 

gl = length(gvalues)
wl = length(wvalues)
N = 4

xticks = [wvalues[1],(N-1)/N,wvalues[end]]
yticks = [round(gvalues[3]/N;digits=2),round((gvalues[1]/N+gvalues[end]/N)/2;digits=2),round(gvalues[end]/N;digits=2)]



replace_infs!(Λ,-1)
pd = plot_phase_diagram(Λ,wvalues,gvalues,N;HEATMAP_DEFAULTS...,xticks=xticks,yticks=yticks)
#annotate!(pd,(N-1)/N,0.5,text(L"w_c", 15,:bottom))


plot!(pd,ns_data[:,2],ns_data[:,1] ./ N,st=:scatter,ms=1,mc=:darkgreen,markerstrokewidth=0.0,markerstrokealpha=0.0,markerstrokecolor=:darkgreen)
plot!(pd,ns_data[:,3],ns_data[:,1] ./ N,st=:scatter,ms=1,mc=:darkgreen,markerstrokewidth=0.0,markerstrokealpha=0.0,markerstrokecolor=:darkgreen)
scatter!([0.7,0.71,0.692],[6.0/N,6.0/N,8.0/N],markershape=:star5,mc=[:gray70, :gray40, :gray10],ms=2,markerstrokecolor=[:gray70, :gray40, :gray10])

lens_xmin = 0.68
lens_xmax = 0.725
lens_ymin = 5.5/N
lens_ymax = 8.2/N

lens_width = lens_xmax - lens_xmin
lens_height = lens_ymax - lens_ymin

lens!([lens_xmin,lens_xmax],[lens_ymin,lens_ymax],inset=(1,bbox(0.6,0.05,0.15,0.55)),
	lw=1.2,
	lc=:gray10,
	clims=(-1/2,1/2),
	framestyle=:box)
	
#pd[2] is the inset
plot!(pd[2],colorbar=false,ticks=false,lw=10)
scatter!(pd[2],[0.7,0.71,0.692],[6.0/N,6.0/N,8.0/N],markershape=:star5,legend=false,mc=[:gray70, :gray40, :gray10],ms=6,markerstrokecolor=[:gray70, :gray40, :gray10])
savefig(pd,output_filename)
pd

