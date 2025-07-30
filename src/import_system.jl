export att_din!,skipsum,create_original_ds,create_ds,generate_initial_condition,rescaled_params

function skipsum(u,i)
    sum = 0.0
    for j in 1:2:length(u)
        sum += u[j]
    end
    return sum - u[i]
end

function att_din!(du,u,p,t)
    n = convert(Int32,p[3])
    for i in 1:2:(2*n)
        du[i] = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(n-1)*skipsum(u,i))
        du[i+1] = u[i] - u[i+1]
    end
end

function att_din(u,p,t)
    n = convert(Int32,p[3])
    du = zeros(length(u))
    for i in 1:2:(2*n)
        du[i] = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(n-1)*skipsum(u,i))
        du[i+1] = u[i] - u[i+1]
    end
    du
end

function att_din_original!(du,u,p,t)
    rp,rc,K,c,α = p
    N = Int(length(u)/2)
    for i in 1:2:(2*N)
        du[i] = rp*u[i]*(1 - rc/K*u[i+1] - c*skipsum(u,i))
        du[i+1] = u[i] - α*u[i+1]
    end
end


diffeq_default =  (abstol = 1.0e-6, reltol = 1.0e-6)
function create_ds(u0,p;diffeq=diffeq_default) CoupledODEs(att_din!,u0,p;diffeq=diffeq) end
function create_integrator(ds,diffeq) integrator(ds;diffeq) end
function create_original_ds(u0,p;diffeq=diffeq_default) CoupledODEs(att_din_original!,u0,p;diffeq=diffeq) end

#convert original params to rescaled
rescaled_params(rp,rc,K,c,α,N) = ( 1/(1 + rc/(α*K*c*(N-1))), rp/α )

#-----------initial conditions-----------------
#define well-controlled initial condition for each system size

function generate_initial_condition(N)
    u0 = zeros(Float64,2N)
    for i in 1:2:2N 
        u0[i] = (i+1)/2
        u0[i+1] = 0.0 
    end
    return u0
end