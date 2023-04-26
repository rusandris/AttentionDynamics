function my_linreg(X,Y)
	x = y = xx = xy = 0
	n=length(X)
	for i in 1:n
		x += X[i]
		xx += X[i]^2
		y += Y[i]
		xy += Y[i]*X[i]
	end
	a=(x*y-n*xy)/(x^2-n*xx)
	b=(x*xy-y*xx)/(x^2-n*xx)

	return (a,b)
end
