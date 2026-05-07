using Revise
#using DynamicalSystems
using AttentionDynamics
#using OrdinaryDiffEq
#using DelimitedFiles
using StatsBase
using LinearAlgebra
using Plots,LaTeXStrings
cd(@__DIR__)

#IDEA: demonstrate fixed point calclulations with figure
#make a plot where M=0, M=1, 1 < M < N and M=N cases are shown 
#plot analytical and numerical eigenvalue spectra for each case 
#discuss also g < N, and g > N cases (separate plots)

function extract_nondegenerate(vals;atol=1e-2)

    eigs_nondeg = [vals[1]]

    for i in 2:length(vals)
        val = vals[i]
        deg = any([isapprox(val, eig,atol=atol) for eig in eigs_nondeg])
        
        if !(deg)
            push!(eigs_nondeg,val)
        end
    end
    return eigs_nondeg
end

function add_NaNs!(vals)
    lmax = maximum(length.(vals))
    for i in 1:length(vals)
        l = length(vals[i])
        push!(vals[i],fill(NaN,lmax-length(vals[i]))...)
    end
end
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])


#------------------------set parameters----------------------------

#N = 4

#N = 5
#g = 5.0
#N = 10
#g = 45.0

#=
w_test = wc1+ (wc2-wc1)/2
diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
ds = create_ds(rand(2N),[g,w_test,N];diffeq=diffeq)
=#

N = 3
ws = [0.0:0.01:1.0;]

gs = [4.0,2.5] 
plots = []
x_axis_setting =:none
y_axis_setting =:auto

ylims = (-5,5)
ylims2 = (-3, 3)
yticks = [ylims[1], 0, ylims[end]]
yticks2 = [ylims2[1], 0, ylims2[end]]
xticks = [0.0,0.5,1.0]

xlims = (0,1)
tickfontsize = 10
guidefontsize = 20
legendfontsize = 13
lw = 3
hlinew = 1
samp = 1 #point sampling
plot_labels = ["a)", "b)"]

critical_w_values = []

for (gi,g) in enumerate(gs)

    wc1 = (N-1)/g
    wc2 = (N-1)/N
    wc3 = (N-1)/(g+N-2)

    push!(critical_w_values,[wc1,wc2,wc3])

    #-------------------------------M = 0---------------------------------

    #eigval spectrum
    #coexistence fixed point 
    u = zeros(2N)
    numerical_spectrums = []
    analytical_spectrums = []

    for (i,w) in enumerate(ws)
        p = [g,w,N]
        @show p

        #build Jacobian
        J = build_jacobian(u,p) 

        #analytical eigenvalues
        analytical_spectrum = fixed_point_eigen_spectrum(N,0,p)
        vals_real = real.(analytical_spectrum)
        vals_imag = imag.(analytical_spectrum)

        push!(analytical_spectrums,vals_real)

    end

    #-----------------------------plot real part------------------------------------
    S_a = transpose(stack(analytical_spectrums))
    #S_n = transpose(stack(numerical_spectrums))

    colors = [:blue :orange]
    labels = [L"\lambda_a" L"\lambda_b"]

    pl0 = plot(xlims=xlims, ylims=ylims,
        xformatter=x_axis_setting,yformatter=y_axis_setting,
        legendcolumns=-1,legendfontsize=legendfontsize,
        guidefontsize=guidefontsize,tickfontsize=tickfontsize,
        yticks=yticks, xticks=xticks)

    #add title to top row figues 
    gi == 1 && plot!(pl0, title=L"^{0}\mathbf{u}^{*}", titlefontsize=16)

    bool_mask = S_a .< 0
    linestyle_mask = similar(bool_mask,Symbol)
    linewidth_mask = similar(bool_mask, Int64)
    for i in 1:length(bool_mask)
        linestyle_mask[i] = bool_mask[i] ? :solid : :solid
        linewidth_mask[i] = bool_mask[i] ? 4 : 1
    end
    #println(linestyle_mask)
    plot!(ws[1:samp:end,1:end],S_a[1:samp:end,1:end],
        ls=linestyle_mask[1:samp:end, 1:end], lw=lw, 
        lc=colors,label = labels)
    #plot!(ws, S_n, st=:scatter, ms=4, mc=colors,
    #    markerstrokewidth=0.0,label = labels_num)
    hline!([0.0],lc=:gray30,ls=:dash,label="",lw=hlinew)

    #add stability annotation 
    annotate!(0.5, 4.5, text("unstable", :gray10, :center, 12))

    annotate!(-0.35, 7.2, text(plot_labels[gi] * " " * L" 1/\epsilon = %$(gs[gi])", :gray10, :left, 16))

    #-------------------------------M = 1---------------------------------

    #eigval spectrum
    #coexistence fixed point 
    u = zeros(2N)
    x_fp(w,M,N) = (N-1)/(N-1 + w*M - w*N) #fixed point value (nonzero)

    numerical_spectrums = []
    analytical_spectrums = []

    for (i,w) in enumerate(ws)
        p = [g,w,N]

        x_fp_value = x_fp(w,1,N) 
        u[[1,N+1]] .= [x_fp_value,x_fp_value] #set both x1 and y1 to x

        #build Jacobian
        J = build_jacobian(u,p)

        #analytical eigenvalues
        analytical_spectrum = fixed_point_eigen_spectrum(N,1,p)
        vals_real = real.(analytical_spectrum)
        vals_imag = imag.(analytical_spectrum)

        push!(analytical_spectrums,vals_real)

    end

    #-----------------------------plot real part------------------------------------
    S_a = transpose(stack(analytical_spectrums))
    #S_n = transpose(stack(numerical_spectrums))

    colors = [:blue :orange :green :green ]
    labels = [L"\lambda_a" L"\lambda_b" L"\lambda_d" ""]
    labels_num = [L"\lambda_a^{num}" L"\lambda_b^{num}" L"\lambda_d^{num}"]

    pl1 = plot(xlims=xlims, ylims=ylims2,
        xformatter=x_axis_setting, yformatter=y_axis_setting,
        legendcolumns=-1, legendfontsize=legendfontsize, legend=:bottomleft,
        guidefontsize=guidefontsize, tickfontsize=tickfontsize,
        yticks=yticks2,xticks=xticks)

    #add title to top row figues 
    gi == 1 && plot!(pl1, title=L"^{1}\mathbf{u}^{*}", titlefontsize=16)

    bool_mask = S_a .< 0
    linestyle_mask = similar(bool_mask,Symbol)
    for i in 1:length(bool_mask)
        linestyle_mask[i] = bool_mask[i] ? :solid : :solid
    end

    plot!(ws[1:samp:end],S_a[1:samp:end,1:end],ls=linestyle_mask[1:samp:end,1:end], lw=lw, lc=colors, label = labels)
    #plot!(ws, S_n, st=:scatter, ms=4, mc=colors,
    #    markerstrokewidth=0.0,label = labels_num)
    hline!([0.0],lc=:gray30,ls=:dash,label="",lw=hlinew)

    #add stability annotation 
    annotate!(0.3, 2.5 , text("unstable", :gray10, :center, 12))
    annotate!(0.85, 2.5 , text("stable", :gray10, :center, 12))
    annotate!(0.85, 1.9, text("dom.", :gray10, :center, 12))

    #add critical value annotation 
    annotate!(critical_w_values[gi][2], -3.6, text(L"w_{\mathrm{d}}", :gray10, :center, 16))

    #-------------------------------1 < M < N ---------------------------------

    #eigval spectrum
    #coexistence fixed point 
    u = zeros(2N)

    numerical_spectrums = []
    analytical_spectrums = []

    for (i,w) in enumerate(ws)
        p = [g,w,N]

        x_fp_value = x_fp(w,2,N) 
        u[[1,N+1]] .= [x_fp_value,x_fp_value] #set both x1 and y1 to fp
        u[[2,N+2]] .= [x_fp_value,x_fp_value] #set both x2 and y2 to fp

        #build Jacobian
        J = build_jacobian(u,p)

        #analytical eigenvalues
        analytical_spectrum = fixed_point_eigen_spectrum(N,2,p)
        vals_real = real.(analytical_spectrum)
        vals_imag = imag.(analytical_spectrum)

        push!(analytical_spectrums,vals_real)

    end

    #-----------------------------plot real part------------------------------------
    add_NaNs!(analytical_spectrums)
    S_a = transpose(stack(analytical_spectrums))
    #S_n = transpose(stack(numerical_spectrums))

    colors = [:blue :orange :green :green :purple :purple]
    labels = [L"\lambda_a" L"\lambda_b" L"\lambda_c" "" L"\lambda_d" ""]
    labels_num = [L"\lambda_a^{num}" L"\lambda_b^{num}" L"\lambda_c^{num}"]

    bool_mask = S_a .< 0
    linestyle_mask = similar(bool_mask,Symbol)
    for i in 1:length(bool_mask)
        linestyle_mask[i] = bool_mask[i] ? :solid : :solid
    end

    @show size(S_a)
    @show size(linestyle_mask)
    pl2 = plot(xlims=xlims, ylims=ylims2,
        xformatter=x_axis_setting,yformatter=y_axis_setting,
        legendcolumns=-1, legendfontsize=legendfontsize, legend=:bottomleft,
        guidefontsize=guidefontsize,tickfontsize=tickfontsize,
        yticks=yticks2,xticks=xticks)

    #add title to top row figues 
    gi == 1 && plot!(pl2, title=L"^{2}\mathbf{u}^{*}", titlefontsize=16)

    plot!(ws[1:samp:end],S_a[1:samp:end,1:end],ls=linestyle_mask[1:samp:end,1:end], lw=lw, lc=colors, label = labels)
    #plot!(ws, S_n, st=:scatter, ms=4, mc=colors,
    #    markerstrokewidth=0.0,label = labels_num)
    hline!([0.0],lc=:gray30,ls=:dash,label="",lw=hlinew)

    #add stability annotation 
    annotate!(0.5, 2.5, text("unstable", :gray10, :center, 12))

    #-------------------------------M = N ---------------------------------

    #eigval spectrum
    #coexistence fixed point 
    u = ones(2N)

    numerical_spectrums = []
    analytical_spectrums = []

    for (i,w) in enumerate(ws)
        p = [g,w,N]

        #build Jacobian
        J = build_jacobian(u,p)

        #analytical eigenvalues
        analytical_spectrum = fixed_point_eigen_spectrum(N,3,p)
        vals_real = real.(analytical_spectrum)
        vals_imag = imag.(analytical_spectrum)

        push!(analytical_spectrums,vals_real)

    end

    #-----------------------------plot real part------------------------------------
    add_NaNs!(analytical_spectrums)
    S_a = transpose(stack(analytical_spectrums))
    #S_n = transpose(stack(numerical_spectrums))

    colors = [:green :green :purple :purple]
    labels = [L"\lambda_c" "" L"\lambda_d" ""]
    labels_num = [L"\lambda_a^{num}" L"\lambda_b^{num}" L"\lambda_c^{num}"]

    bool_mask = S_a .< 0
    linestyle_mask = similar(bool_mask,Symbol)
    for i in 1:length(bool_mask)
        linestyle_mask[i] = bool_mask[i] ? :solid : :solid
    end

    pl3 = plot(xlims=xlims, ylims=ylims2,
        xformatter=x_axis_setting,yformatter=y_axis_setting,
        legendcolumns=-1,legendfontsize=legendfontsize,
        guidefontsize=guidefontsize,tickfontsize=tickfontsize,
        yticks=yticks2,xticks=xticks,legend=:bottomleft)

    #add title to top row figues 
    gi == 1 && plot!(pl3, title=L"^{N}\mathbf{u}^{*}", titlefontsize=16)

    plot!(ws[1:samp:end],S_a[1:samp:end,1:end],ls=linestyle_mask[1:samp:end,1:end], lw=lw, lc=colors, label = labels)
    #plot!(ws, S_n, st=:scatter, ms=4, mc=colors,
    #    markerstrokewidth=0.0,label = labels_num)
    hline!([0.0],lc=:gray30,ls=:dash,label="",lw=hlinew)

    #add stability annotation 
    annotate!(0.3, 2.5, text("stable", :gray10, :center, 12))
    annotate!(0.3, 1.9, text("coex.", :gray10, :center, 12))

    if gi == 1
        annotate!(0.73, 2.5, text("unstable", :gray10, :center, 12))
    elseif gi == 2
        annotate!(0.87, 2.5, text("unstable", :gray10, :center, 12))
    end

    #add critical value annotation 
    annotate!(critical_w_values[gi][1], -3.6, text(L"w_{\mathrm{c}}", :gray10, :center, 16))
    annotate!(critical_w_values[gi][2], -3.6, text(L"w_{\mathrm{d}}", :gray10, :center, 16))

    #pl_g1 = plot(pl0,pl1,pl2,pl3,layout=(1,4),size=(1200,300))
    push!(plots,[pl0,pl1,pl2,pl3])

end

#edit figures
for i in 1:length(plots)
    pls = plots[i]
    ws_crit = critical_w_values[i]

    for j in 1:length(pls)

        #add vertical lines
        if j == 2 
            vline!(pls[j],[ws_crit[2]],lc=:gray10,lw=2,ls=:dot,label="") #label=L"w^{(2)}_c")
        elseif j == 4
            vline!(pls[j],[ws_crit[1]],lc=:gray10,lw=2,ls=:dash,label="") #,label=L"w^{(1)}_c")
            vline!(pls[j], [ws_crit[2]], lc=:gray10, lw=2, ls=:dot, label="") #label=L"w^{(2)}_c")

            if i == 2
                vline!(pls[j], [ws_crit[2]], lc=:gray10, lw=2, ls=:dot, label="") #label=L"w^{(2)}_c")
                #vline!(pls[j],[ws_crit[1]],lc=:gray10,lw=2,ls=:dash,label="") #,label=L"w^{(1)}_c")
                #vline!(pls[j],[ws_crit[3]],lc=:gray50,lw=2,ls=:dash,label="") #,label=L"w^{(1)}_c")
            end
        end
        
        #set y-axis attributes
        if j == 1
            plot!(pls[j],yformatter=:auto,ylabel=L"Re(\lambda)")
        end

        #set x-axis attributes
        if i == 2
            plot!(pls[j],xformatter=:auto,xlabel=L"w")
        end

        #draw shading
        
        #draw shading for M=1 cases
        if j == 2
            width = ws[end]-ws_crit[2]
            height = ylims[end] - ylims[1]
            plot!(rectangle(width,height,ws_crit[2],ylims[1]), 
                opacity=.2,color=:gray50,label="")
        end
            
        #draw shading for M=N cases

        if j == 4
            width = ws_crit[i] - ws[1]
            height = ylims[end] - ylims[1]
            plot!(rectangle(width,height,0,ylims[1]), 
                opacity=.2,color=:gray50,label="")
        end


        #add margins to first row figs 
        if i == 1
            plot!(pls[j], top_margin=8Plots.mm)
        end

        #add margins to second row figs 
        if i == 2
            plot!(pls[j], top_margin=10Plots.mm, bottom_margin=5Plots.mm)
        end

        #add left margins to first column figs 
        if j == 1
            plot!(pls[j], left_margin=8Plots.mm, top_margin=3Plots.mm)
        end

    end
end

pl_final = plot(plots[1]...,plots[2]...,layout = (2,4),
    size=(1200, 500), dpi=400)
    #left_margin=-1Plots.mm,dpi=300,
    #right_margin=1Plots.mm,bottom_margin=-1Plots.mm)

figs_dir = "../../paper_figs/"
savefig(pl_final,figs_dir*"eigen_spectrum_fig_N3_g2.5_4.png")
display(pl_final)














#=

#dominance fixed point: doesn't matter which one we choose 
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=1
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)

u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,1,p)

#M=2 fixed point: check one of those
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=2
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
u[[2,N+2]] .= [s,s] #set both x2 and y2 to s

J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,2,p)

#M=3 fixed point: check one of those
u = zeros(2N)
p = [g,0.7,N]
w = p[2]
M=3
w_M = (w*M)/(N-1)
w_N = (w*N)/(N-1)
s = 1/(1 + w_M - w_N) #fixed point value (nonzero)
u[[1,N+1]] .= [s,s] #set both x1 and y1 to s
u[[2,N+2]] .= [s,s] #set both x2 and y2 to s
u[[3,N+3]] .= [s,s] #set both x3 and y3 to s

J = build_jacobian(u,p)
#test if results are correct
numerical_spectrum =  eigvals!(J)
analytical_spectrum = fixed_point_eigen_spectrum(N,3,p)

#parameter dependence -> look for bifurcations of limit cycles in terms of w

ws = [0.3:0.01:0.9;]

#reals,imags=eigen_spectrum(J,ws)
λ_specs_coex = []
λ_specs_dom = []
λ_specs_M2 = []
λ_specs_M3 = []
λ_specs_M4 = []
lyaps = []
max_crs = [] #contraction rate maxima

for (i,w) in enumerate(ws)
    p = [g,w,N]
    spec_coex = fixed_point_eigen_spectrum(N,N,p) #coexistence spectrum
    spec_dom = fixed_point_eigen_spectrum(N,1,p) #dominance spectrum
    spec_M2 = fixed_point_eigen_spectrum(N,2,p) #M=2 spectrum
    spec_M3 = fixed_point_eigen_spectrum(N,3,p) #M=3 spectrum
    spec_M4 = fixed_point_eigen_spectrum(N,4,p) #M=4 spectrum

    push!(λ_specs_coex,spec_coex)
    push!(λ_specs_dom,spec_dom)
    push!(λ_specs_M2,spec_M2)
    push!(λ_specs_M3,spec_M3)
    push!(λ_specs_M4,spec_M4)
    
    
    if w > wc1
        set_parameters!(ds,p)
        push!(lyaps,lyapunov(ds,1000;Ttr=100))
        traj,ts = trajectory(ds,200;Ttr=10000,Δt=1e-2)
        crs = [contraction_rate(traj[i],p) for i in 1:length(traj)]
        min_cr,max_cr = extrema(crs)
        @show max_cr
        push!(max_crs,max_cr)
    end
    
    @show w
end

#extract real and imaginary parts 
realparts_coex = stack(real.(λ_specs_coex))'
imagparts_coex = stack(imag.(λ_specs_coex))'
realparts_dom = stack(real.(λ_specs_dom))'
imagparts_dom = stack(imag.(λ_specs_dom))'
realparts_M2 = stack(real.(λ_specs_M2))'
imagparts_M2 = stack(imag.(λ_specs_M2))'
realparts_M3 = stack(real.(λ_specs_M3))'
imagparts_M3 = stack(imag.(λ_specs_M3))'
realparts_M4 = stack(real.(λ_specs_M4))'
imagparts_M4 = stack(imag.(λ_specs_M4))'


ylims_spec = [-6.0,3.0]

#plot coexistence fixed point spectrum
pl = plot(ws,realparts_coex[:,[1,2]],ylabel="eigen spectrum",xlabel=L"w",legend=:outertopright,
    title="N=$(N),g=$(g)",ylims=ylims_spec,lc=:orange,lw=2,label="real type 3",guidefontsize=16,dpi=300)
plot!(ws,realparts_coex[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,label="real type 4")
plot!(ws,imagparts_coex[:,[1,2]],ylims=ylims_spec,lc=:orange,lw=2,ls=:dash,label="imag type 3")
plot!(ws,imagparts_coex[:,[3,4]],ylims=ylims_spec,lc=:green,lw=2,ls=:dash,label="imag type 4")

#plot dominance fixed point spectrum 
#plot!(ws,realparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 dom")
#plot!(ws,imagparts_dom[:,[3,4]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="real type 4 dom")

#plot M=2 fixed point spectrum
plot!(ws,realparts_M2[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,label="real type 3 M=2")
plot!(ws,imagparts_M2[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,ls=:dash,label="imag type 3 M=2")
plot!(ws,realparts_M2[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 M=2")
plot!(ws,imagparts_M2[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="imag type 4 M=2")

#plot M=3 fixed point spectrum
plot!(ws,realparts_M3[:,[3,4]],ylims=ylims_spec,lc=:purple,lw=2,label="real type 3 M=3")
plot!(ws,imagparts_M3[:,[3,4]],ylims=ylims_spec,lc=:purple,lw=2,ls=:dash,label="imag type 3 M=3")
plot!(ws,realparts_M3[:,[5,6]],ylims=ylims_spec,lc=:brown,lw=2,label="real type 4 M=3")
plot!(ws,imagparts_M3[:,[5,6]],ylims=ylims_spec,lc=:brown,lw=2,ls=:dash,label="imag type 4 M=3")

#plot M=4 fixed point spectrum
plot!(ws,realparts_M4[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,label="real type 3 M=4")
plot!(ws,imagparts_M4[:,[3,4]],ylims=ylims_spec,lc=:red,lw=2,ls=:dash,label="imag type 3 M=4")
plot!(ws,realparts_M4[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,label="real type 4 M=4")
plot!(ws,imagparts_M4[:,[5,6]],ylims=ylims_spec,lc=:blue,lw=2,ls=:dash,label="imag type 4 M=4")



vline!([wc1],lc=:gray10,lw=2,ls=:dash,label=L"w^{(1)}_c")
vline!([wc2],lc=:gray10,lw=2,ls=:dash,label=L"w^{(2)}_c")
#vline!([0.55],lc=:gray10,lw=2,ls=:dash,label=L"0.55")
zeropoint_M2 = findlast(x -> x != 0.0, imagparts_M2[:,3])
zeropoint_coex = findlast(x -> x != 0.0, imagparts_coex[:,1])

#ws[zeropoint_M2]
#ws[zeropoint_coex]
#vline!([ws[zeropoint_M2]],lc=:green,lw=2,ls=:dash,label=L"pd")
#vline!([ws[zeropoint_coex]],lc=:green,lw=2,ls=:dash,label=L"pd")

plot!(ws[25:end],lyaps,lc=:gray10,lw=2,label=L"\lambda")
plot!(ws[25:end],max_crs,lc=:purple,lw=2,label="cr max")

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(ws[end]-ws[1],6,ws[1],-6), opacity=.2,fillc=:gray50,label="")

savefig(pl,"../../eigen_spectrum_coexistence_fp_N$(N).png")
=#