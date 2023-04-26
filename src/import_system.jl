function skipsum(a,i)
	x = a[1:2:end]
    return sum(x) - a[i]
end

function att_din!(du,u,p,t)
    n = convert(Int32,p[3])
    for i in 1:2:(2*n)
        du[i] = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(n-1)*skipsum(u,i))
        du[i+1] = u[i] - u[i+1]
    end
end

function create_ds(u0,p) ContinuousDynamicalSystem(att_din!,u0,p) end
function create_integrator(ds,diffeq) integrator(ds;diffeq) end



