export build_jacobian, build_jacobian_symLV, eigval_maxrealparts, eigen_spectrum
export fixed_point_eigen_spectrum,contraction_rate

function skipsum(a,i,N)
	x = @view a[1:N]
    return sum(x) - x[i]
end

#build jacobian for attdin
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

#build jacobian for symmetric classic LV
function build_jacobian_symLV(u,p,M)
    N = length(u)
    r, w1, w2 = p
    J = zeros(N, N)

    for i in 1:N

        #elements of J
        Ai = A(u, p, i)
        Bi = B(u, p, i)

        for j in 1:N            
            if i == j 
                J[i,j] = Ai
            else 
                J[i, j] = Bi
            end
        end
    end
    return J
end

function A(u,p,i)
    r, w1, w2 = p
    return r*(1 - w2*skipsum(u,i;step=1) - 2*w1*u[i])
end

function B(u,p,i)
    r, w1, w2 = p
    return -r*u[i]*w2
end

#analytical symmetric LV eigenspectrum
function fixed_point_eigen_spectrum(N,M,p)
    r, w1, w2 = p

    x = 1/(w1+(M-1)*w2)
    A = r*(1 - ( (M-1)*w2 + 2w1)*x )
    A0 = r * (1 - M * w2 * x)
    B = -r*w2*x

    λ1 = A + B*(M-1)
    λ2 = A - B 
    λ3 = A0
    
    return [λ1,fill(λ2,M-1)...,fill(λ3,N-M)...]

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
	g,w,N = p
	
	#need expressions, convenient notations
	w_M = (w*M)/(N-1)
	w_N = (w*N)/(N-1)
	s = 1/(1 + w_M - w_N) #fixed point value (nonzero)

	a2 = -(g*w*s)/(N-1)
	b = g*(w-1)*s
	a3 = -(M-1)*a2


	eigen_spectrum = []

	#eigenvalues 
	λ1 = -1.0 #type 1
	λ2 = g*(1 - w_M*s) #type 2
	#type 3
	λ3_1 = 0.5*( -(1 + a2) + sqrt(Complex((1+a2)^2 + 4(b-a2))) )
	λ3_2 = 0.5*( -(1 + a2) - sqrt(Complex((1+a2)^2 + 4(b-a2))) )
	# type 4
	λ4_1 = 0.5*( -(1 + a3) + sqrt(Complex((1+a3)^2 + 4(b-a3))) )
	λ4_2 = 0.5*( -(1 + a3) - sqrt(Complex((1+a3)^2 + 4(b-a3))) )

	#origin case 
	if M == 0
		push!(eigen_spectrum,λ1,λ2)
		return eigen_spectrum
	end

	#dominance case
	if M == 1
		push!(eigen_spectrum,λ1,λ2,λ4_1,λ4_2)
		return eigen_spectrum
	end

	if M == N 
		push!(eigen_spectrum,λ3_1,λ3_2,λ4_1,λ4_2)
		return eigen_spectrum
	else #all other cases
		push!(eigen_spectrum,λ1,λ2,λ3_1,λ3_2,λ4_1,λ4_2)
	end
	return eigen_spectrum
end

#phase space volume contraction rate
function contraction_rate(u,p)
	g,w,N = p
	sum_x = sum(u[1:2:end])
	sum_y = sum(u[2:2:end])
	return g*N - g* sum_y - g*w*(sum_x-sum_y) - N
end
