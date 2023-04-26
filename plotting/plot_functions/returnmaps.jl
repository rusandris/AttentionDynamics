function plot_returnmap(ts;filename,ms,denseticks=false,xlabel,ylabel)

	if denseticks
		rm = plot(ts[1:end-1],ts[2:end],
		st=:scatter,
		mc=:gray10,
		xlabel=xlabel,
		ylabel=ylabel,
		legend=false,
		xguidefontsize=18,
		yguidefontsize=18,
		tickfontsize=10,
		ms=ms,
		framestyle = :box,
		ticks = 0:maximum(ts),
		dpi = 300)
	else
		rm = plot(ts[1:end-1],ts[2:end],
		st=:scatter,
		mc=:gray10,
		xlabel=xlabel,
		ylabel=ylabel,
		legend=false,
		xguidefontsize=18,
		yguidefontsize=18,
		tickfontsize=10,
		ms=ms,
		framestyle = :box,
		dpi = 300)
	
	end
	
	savefig(rm,filename)
end

function plot_color_returnmap(ts;xlabel,ylabel,ms,color=:thermal,denseticks=false,filename)
	transitions = collect(eachrow(hcat(ts[1:end-1],ts[2:end])))
	cnts = counter(transitions)
	
	if denseticks == true
		pl = plot(ticks = 0:maximum(ts))
	else
		pl = plot()
	end
	
	for trn in transitions
		prob = cnts[trn] ./ length(transitions)
		plot!(pl,[trn[1]],[trn[2]],
		st=:scatter,
		ms=ms,
		zcolor = prob,
		color=color,
		xlabel=xlabel,
		ylabel=ylabel,
		legend=false,
		dpi=300,
		framestyle = :box,
		xguidefontsize=16,
		yguidefontsize=16,
		markerstrokewidth=0) 
	end
	
	savefig(pl,filename)
end

function color_returnmap(ts;xlabel,ylabel,ms,color=:thermal,mode=:lin,denseticks=false,filename)
	states = @view ts[1:end-1]
	next_states = @view ts[2:end]
	
	C = spzeros(1023,1023)
	
	for i in 1:length(states)
		C[states[i],next_states[i]] += 1
	end
	
	
	
	if denseticks == true
		pl = plot(ticks = 0:maximum(ts))
	else
		pl = plot()
	end
	
	if mode ==:lin
		for i in 1:length(states)
			prob = C[states[i],next_states[i]] / sum(C)
			plot!(pl,[states[i]],[next_states[i]],
			st=:scatter,
			ms=ms,
			zcolor = prob,
			color=color,
			xlabel=xlabel,
			ylabel=ylabel,
			legend=false,
			dpi=300,
			markerstrokewidth=0,
			framestyle = :box,
			xguidefontsize=16,
			yguidefontsize=16) 
		end
	elseif mode ==:log
		for i in 1:length(states)
			prob = C[states[i],next_states[i]] / sum(C)
			plot!(pl,[log2(states[i])],[log2(next_states[i])],
			st=:scatter,
			ms=ms,
			zcolor = prob,
			color=color,
			xlabel=xlabel,
			ylabel=ylabel,
			legend=false,
			dpi=300,
			markerstrokewidth=0,
			framestyle = :box,
			xguidefontsize=16,
			yguidefontsize=16) 
		end
	end
	
	savefig(pl,filename)
	
end

function returnmap_heatmap(ts;xlabel,ylabel,color,filename,rev=false)
	states = @view ts[1:end-1]
	next_states = @view ts[2:end]
	
	C = spzeros(1023,1023)
	
	for i in 1:length(states)
		C[states[i],next_states[i]] += 1
	end
	
	
	pl = plot(Matrix(C),st=:heatmap,
	framestyle = :box,
	xguidefontsize=16,
	yguidefontsize=16,
	color=color,
	rev=false,
	xlabel=xlabel,
	ylabel=ylabel,
	dpi = 300)
	
	savefig(pl,filename)
end

function number_of_flips(state,next_state)
	nr_flips = 0 
	for i in 1:length(state)
		if state[i] != next_state[i]
			nr_flips += 1
		end
	end
	nr_flips
end


function returnmap_flips(binary_ts,integer_ts,colors;ms,xlabel,ylabel,return_distribution,filename1="maptest.png",filename2="flipdistribtest.png")
	bin_states = @view binary_ts[2:end-1]
	bin_next_states = @view binary_ts[3:end]
	int_states = @view integer_ts[2:end-1]
	int_next_states = @view integer_ts[3:end]
	
	flips = []
	
	pl = plot() 
	
	for i in 1:length(bin_states)
		nr_flips = number_of_flips(bin_states[i],bin_next_states[i])
		
		if return_distribution
			push!(flips,nr_flips)
		end
		
		
		label = "$nr_flips"
		
		labels = [x[:label] for x in pl.series_list]
		if label in labels
			label = ""
		end
		
		plot!(pl,[int_states[i]],[int_next_states[i]],
			st=:scatter,
			ms=ms,
			color=colors[nr_flips+1],
			xlabel=xlabel,
			ylabel=ylabel,
			label=label,
			dpi=300,
			markerstrokewidth=0,
			framestyle = :box,
			xguidefontsize=16,
			yguidefontsize=16) 
	end
	
	if return_distribution
	 	bins = [minimum(flips):maximum(flips)+1;] .- 0.5
	 	
		h = histogram(flips,bins=bins,
		dpi=300,
		normalize=:probability,
		legend=false,
		xticks = minimum(flips):maximum(flips),
		ylabel="P(number of flips)",
		xlabel="number of flips",
		tickfontsize=10,
		xguidefontsize=18,
		yguidefontsize=18)
		
		savefig(pl,filename1)
		savefig(h,filename2)
	else
		savefig(pl,filename1)
	end
end



function plot_recurrence(ts;filename,atol=0,ms=0.2)
	rm = isapprox.(reverse(ts), ts',atol=atol)
	
	rp = spy(rm,
	color = cgrad(:greys),
	legend=:none,
	ticks=:none,
	framestyle=:box,
	ms = ms,
	dpi = 300)
	 
	savefig(rp,filename)
	
end

