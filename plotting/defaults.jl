using Plots

guide_fontsize = 18
tick_fontsize = 12
legend_fontsize = 10

DEFAULTS = (xguidefontsize=guide_fontsize,yguidefontsize=guide_fontsize,legend=false,legendfontsize=legend_fontsize,tickfontsize=tick_fontsize,aspect_ratio=:auto,lw=1,la=1,ls=:solid,st=:path,mc=:gray10,ms=2,grid=false,widen=true,framestyle=:box,right_margin = 15Plots.mm,yguidefontrotation=0,dpi=300)

HEATMAP_DEFAULTS = (xguidefontsize=guide_fontsize,yguidefontsize=guide_fontsize,
    c=cgrad(:seismic, rev = false), clims=(-1/2,1/2),colorbar=true,framestyle=:box,
    lw=4,lc=:gray10,xlabel=L"$w$",ylabel=L"$g/(N-1)$",colorbar_title=L"$\lambda$",xticks=0:0.2:1,
    legendfontsize=legend_fontsize,legend=:bottomleft,colorbar_titlefontsize=guide_fontsize,tickfontsize=tick_fontsize,
    dpi=300,right_margin = 15Plots.mm)

OD_DEFAULTS =  (xlabel=L"$w$",ylabel = L"$x_1$",st=:scatter,ms=0.1,xguidefontsize=guide_fontsize,yguidefontsize=guide_fontsize,tickfontsize=tick_fontsize,mc=:gray10,widen=true,legend=false,dpi=300)

THREED_DEFAULTS = (xguidefontsize=guide_fontsize/2,yguidefontsize=guide_fontsize/2,zguidefontsize=guide_fontsize/2,tickfontsize=tick_fontsize/2,lw=5,legend=false,dpi=300)

BIF_DEFAULTS = (xguidefontsize=guide_fontsize,yguidefontsize=guide_fontsize,legend=false,legendfontsize=legend_fontsize,
tickfontsize=tick_fontsize,grid=false,framestyle=:box,lw=3,ms=8,lc=:gray10,right_margin = 15Plots.mm,yguidefontrotation=0,dpi=300)

BIF_COLORS = Dict(:pd => :orange,:ns => :green, :hopf => :red,:bp => :gray50,:nd => :purple,:d => :gray50,:dhopf => :red)
BIF_NAMES = Dict(:pd => "PD",:ns => "NS", :hopf => "H",:bp => "BP",:nd => "?",:d => "D",:dhopf => "DH")

