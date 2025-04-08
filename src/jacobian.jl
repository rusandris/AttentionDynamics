export build_jacobian,eigval_maxrealparts,eigen_spectrum

function skipsum(a,i,N)
	x = @view a[1:N]
    return sum(x) - x[i]
end

function build_jacobian(u,p)
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

#calc max real part of eigvals numerically
function eigval_maxrealparts(fixed_points,w_values,p)
	l = length(fixed_points[1])
	allmax_realparts = []
	
	for fp in fixed_points
		max_realparts = []
		for w in w_values
			p[2] = w 
			J = build_jacobian(fp,p)
			eigen_values = eigvals!(J)
			push!(max_realparts,maximum(real.(eigen_values)))
		end
		push!(allmax_realparts,max_realparts)
	end
	
	allmax_realparts
end

#calc eigvals of J numerically
function eigen_spectrum(J,w_values)
	realparts = []
	imagparts = []
	
	for w in w_values
		eigen_values = eigvals!(J)
		push!(realparts,real.(eigen_values))
		push!(imagparts,imag.(eigen_values))
	end
	return realparts,imagparts
end

#TODO: implement formulas for analytic eigvals
#extract real parts and imaginary parts

function fixed_point_eigen_spectrum(N,M,p)


end

