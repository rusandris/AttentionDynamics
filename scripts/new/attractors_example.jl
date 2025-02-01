using DynamicalSystems
using Plots,LaTeXStrings

function attention_dynamics_rule(du,u,p,t)
    g,w = p
	N = Int(length(u)/2)

    @inbounds for i in 1:N
    	sum_ui = sum(u[1:N]) - u[i]
        du[i] = g*u[i]*(1 - (1-w)*u[i+N] - w/(N-1) *  sum_ui)
        du[i+N] = u[i] - u[i+N]
    end

end

w = 0.7; g = 18
p = [g,w]
u0 = rand(20)

#create regular ds object (6D dynamical system)
ds = CoupledODEs(attention_dynamics_rule,u0,p)
pmap = PoincareMap(ds,(2,10.0))
#=

#-------------------try attractor finding---------------------------
grid = (
    range(0.0, 5.0; length = 150), # x1
    range(0.0, 5.0; length = 150), # x2
    range(0.0, 5.0; length = 150), # x3
    range(0.0, 3.0; length = 150), # y1
    range(0.0, 3.0; length = 150), # y2
    range(0.0, 3.0; length = 150), # y3
)

mapper = AttractorsViaRecurrences(ds, grid;
    consecutive_recurrences = 1000, attractor_locate_steps = 2000,
    consecutive_lost_steps = 100,
)

attractors = extract_attractors(mapper)
=#