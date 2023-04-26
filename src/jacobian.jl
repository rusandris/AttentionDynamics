function skipsum(a,i,N)
	x = @view a[1:N]
    return sum(x) - x[i]
end

function jacobian(u,p)
	g,w,N = p
	N = Int(N)
	
	XX = zeros(N,N)
	XY = zeros(N,N)
	YX = zeros(N,N)
	YY = zeros(N,N)
	
	
	for i in 1:N
		for j in 1:N
			if i == j
				XX[i,j] = g*(1 - (1 - w)*u[i+N] - w/(N-1)*skipsum(u,i,N))
				XY[i,j] = -g*(1-w)*u[i]
				
				YX[i,j] = 1
				YY[i,j] = -1 
			else
				XX[i,j] = -g*w/(N-1)*u[i]
				XY[i,j] = 0
				
				YX[i,j] = 0
				YY[i,j] = 0
			end
		end
	end
	
	J = vcat(hcat(XX,XY),hcat(YX,YY))
	
end

function eigval_curves(fixed_points,w_values,p)
	l = length(fixed_points[1])
	allmax_realparts = []
	
	for fp in fixed_points
		max_realparts = []
		for w in w_values
			p[2] = w 
			J = jacobian(fp,p)
			eigen_values = eigvals!(J)
			push!(max_realparts,maximum(real.(eigen_values)))
		end
		push!(allmax_realparts,max_realparts)
	end
	
	allmax_realparts
end


