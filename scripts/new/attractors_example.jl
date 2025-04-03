using Revise
using DynamicalSystems
using AttentionDynamics
using OrdinaryDiffEq
using DelimitedFiles
using Plots,LaTeXStrings

diffeq = (alg = Tsit5(),dt = 1e-3,adaptive=true,dtmax=1e-2,maxiters=1e15)
N = 3
g = 5
wc2 = (N-1)/N
wc1 = (N-1)/g
w = 0.5
ds = create_ds(rand(2N),[g,w,N];diffeq=diffeq)

#-------------------try attractor finding---------------------------
grid = (
    range(0.0, 5.0; length = 3), # x1
    range(0.0, 5.0; length = 3), # x2
    range(0.0, 5.0; length = 3), # x3
    range(0.0, 3.0; length = 3), # y1
    range(0.0, 3.0; length = 3), # y2
    range(0.0, 3.0; length = 3), # y3
)

mapper = AttractorsViaRecurrences(ds, grid;Ttr=5000,
    consecutive_recurrences = 5000, attractor_locate_steps = 10000,
    consecutive_lost_steps = 100,
)

#attractors = extract_attractors(mapper)

