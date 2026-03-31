export att_din!, att_din_original!, skipsum, create_original_ds, att_din_original_noisy!, create_original_mod_ds
export create_ds, generate_initial_condition, rescaled_params, rescaling_factors
export create_ds_noisy, create_ds_prop_noisy, create_ds_alle, create_ds_symLV, create_original_ds_noisy,
    create_original_varied_delay, create_original_varied_rate, 
    create_original_varied_rate_delay, create_original_varied_params,
    create_original_asymmetric_interaction

#sum skipping the i-th element
function skipsum(u,i;step=2)
    sum = 0.0
    for j in 1:step:length(u)
        sum += u[j]
    end
    return sum - u[i]
end

function gensum(u, C, idx; step=2)
    return sum(C[idx,:] .* u[1:step:end])
end

#system definition (in-place)
function att_din!(du,u,p,t)
    n = convert(Int32,p[3])
    for i in 1:2:(2*n)
        du[i] = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(n-1)*skipsum(u,i))
        du[i+1] = u[i] - u[i+1]
    end
end

#system definition (out-of-place)
function att_din(u,p,t)
    n = convert(Int32,p[3])
    du = zeros(length(u))
    for i in 1:2:(2*n)
        du[i] = p[1]*u[i]*(1 - (1 - p[2])*u[i+1] - p[2]/(n-1)*skipsum(u,i))
        du[i+1] = u[i] - u[i+1]
    end
    du
end

#system definition (in-place) with additional uniform noise
function att_din_noise!(du, u, p, t)
    n = convert(Int32, p[3])
    for i in 1:2:(2*n)
        du[i] = p[1] * u[i] * (1 - (1 - p[2]) * u[i+1] - p[2] / (n - 1) * skipsum(u, i)) + p[4] * rand()
        du[i+1] = u[i] - u[i+1]
    end
end

#system definition (in-place) with additional uniform noise proportional to x
function att_din_prop_noise!(du, u, p, t)
    n = convert(Int32, p[3])
    for i in 1:2:(2*n)
        du[i] = p[1] * u[i] * (1 - (1 - p[2]) * u[i+1] - p[2] / (n - 1) * skipsum(u, i)) + p[4] * rand() * u[i]
        du[i+1] = u[i] - u[i+1]
    end
end

#system definition (in-place) with additional uniform noise
function att_din_alle!(du, u, p, t)
    n = convert(Int32, p[3])
    for i in 1:2:(2*n)
        du[i] = p[1] * u[i] * (1 - (1 - p[2]) * u[i+1] - p[2] / (n - 1) * skipsum(u, i)) * (u[i] / p[4] - 1)
        du[i+1] = u[i] - u[i+1]
    end
end

#system definition no delay (classic symmetric LV)
function symmetricLV!(du, u, p, t)
    n = length(u)
    r, w1, w2 = p
    for i in 1:n
        du[i] = r * u[i] * (1 - w1 * u[i] - w2 * skipsum(u, i; step=1))
    end
end

#system definition (in-place) with original parameters
function att_din_original!(du, u, p, t)
    rp, rc, K, c, α = p
    N = Int(length(u) / 2)
    for i in 1:2:(2*N)
        du[i] = rp * u[i] * (1 - rc / K * u[i+1] - c * skipsum(u, i))
        du[i+1] = u[i] - α * u[i+1]
    end
end

#system definition (in-place) with varied delay
function att_din_original_varied_delay!(du, u, p, t)
    rp, rc, K, c, αs = p

    #αs vector of delay exponents

    N = Int(length(u) / 2)
    for (idx,i) in enumerate(1:2:(2*N))
        du[i] = rp * u[i] * (1 - rc / K * u[i+1] - c * skipsum(u, i))
        du[i+1] = u[i] - αs[idx] * u[i+1]
    end
end

#system definition (in-place) with varied prod rates
function att_din_original_varied_rate!(du, u, p, t)
    rps, rc, K, c, α = p

    #rps vector of prod rates

    N = Int(length(u) / 2)
    for (idx, i) in enumerate(1:2:(2*N))
        du[i] = rps[idx] * u[i] * (1 - rc / K * u[i+1] - c * skipsum(u, i))
        du[i+1] = u[i] - α * u[i+1]
    end
end

#system definition (in-place) with varied prod rates and delay
function att_din_original_varied_rate_delay!(du, u, p, t)
    rps, rc, K, c, αs = p

    #αs vector of delay exponents
    #rps vector of prod rates

    N = Int(length(u) / 2)
    for (idx, i) in enumerate(1:2:(2*N))
        du[i] = rps[idx] * u[i] * (1 - rc / K * u[i+1] - c * skipsum(u, i))
        du[i+1] = u[i] - αs[idx] * u[i+1]
    end
end

#system definition (in-place) with varied prod rates and delay
function att_din_original_varied_params!(du, u, p, t)
    rps, rc, K, cs, αs = p

    #αs vector of delay exponents
    #rps vector of prod rates
    #cs two values: high and low

    N = Int(length(u) / 2)
    for (idx, i) in enumerate(1:2:(2*N)) 

        du[i] = rps[idx] * u[i] * (1 - rc / K * u[i+1] - cs[idx] * skipsum(u, i))
        du[i+1] = u[i] - αs[idx] * u[i+1]
    end
end


#system definition (in-place) with varied prod rates and delay
function att_din_asymmetric_interaction!(du, u, p, t)
    rp, rc, K, C, α = p

    N = Int(length(u) / 2)
    for (idx,i) in enumerate(1:2:(2*N))
        du[i] = rp * u[i] * (1 - rc / K * u[i+1] - gensum(u,C,idx))
        du[i+1] = u[i] - α * u[i+1]
    end

end

#system definition (in-place) with original parameters with added noise
function att_din_original_noisy!(du, u, p, t)
    rp, rc, K, c, α, σ = p
    N = Int(length(u) / 2)
    for i in 1:2:(2*N)
        du[i] = rp * u[i] * (1 - rc / K * u[i+1] - c * skipsum(u, i)) + σ * rand()
        du[i+1] = u[i] - α * u[i+1]
    end
end

#system definition (in-place) with original parameters modified carrying capacity
function att_din_original_mod!(du, u, p, t)
    rp, rc, K, c, α = p
    N = Int(length(u) / 2)
    for i in 1:2:(2*N)
        du[i] = rp * u[i] * (K - rc * u[i+1] - c * skipsum(u, i))
        du[i+1] = u[i] - α * u[i+1]
    end
end

diffeq_default =  (abstol = 1.0e-6, reltol = 1.0e-6)
function create_ds(u0,p;diffeq=diffeq_default) CoupledODEs(att_din!,u0,p;diffeq=diffeq) end
function create_integrator(ds,diffeq) integrator(ds;diffeq) end
function create_original_ds(u0,p;diffeq=diffeq_default) CoupledODEs(att_din_original!,u0,p;diffeq=diffeq) end
function create_original_ds_noisy(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_noisy!, u0, p; diffeq=diffeq) end
function create_original_varied_delay(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_varied_delay!, u0, p; diffeq=diffeq) end
function create_original_varied_rate(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_varied_rate!, u0, p; diffeq=diffeq) end
function create_original_varied_rate_delay(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_varied_rate_delay!, u0, p; diffeq=diffeq) end
function create_original_varied_params(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_varied_params!, u0, p; diffeq=diffeq) end
function create_original_asymmetric_interaction(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_asymmetric_interaction!, u0, p; diffeq=diffeq) end
function create_original_mod_ds(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_original_mod!, u0, p; diffeq=diffeq) end
function create_ds_noisy(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_noise!, u0, p; diffeq=diffeq) end
function create_ds_prop_noisy(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_prop_noise!, u0, p; diffeq=diffeq) end
function create_ds_alle(u0, p; diffeq=diffeq_default) CoupledODEs(att_din_alle!, u0, p; diffeq=diffeq) end
function create_ds_symLV(u0, p; diffeq=diffeq_default) CoupledODEs(symmetricLV!, u0, p; diffeq=diffeq) end


#convert original params to rescaled
rescaled_params(rp,rc,K,c,α,N) = ( 1/(1 + rc/(α*K*c*(N-1))), rp/α )
rescaling_factors(rp, rc, K, c, α, N) = ( α*K/(rc+α*K*c*(N-1)), K/(rc + α*K*c*(N-1)) )
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