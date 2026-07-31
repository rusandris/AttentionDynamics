# AttentionDynamics.jl
Source code for dynamics of collective attention based on the paper Lorenz-Spreen et al. (2019).
Functionality includes:
* numerical integration of the high-dimensional Lotka-Volterra system in various forms
* evaluation of Lyapunov exponents
* charting statistics: definition of popularity lifetimes from trajectories
* plotting scripts for figures

WIP, stay tuned for fully organised, documented version

## Installation
Download and install automatically to your local/global `julia` environment with:
```julia
] add https://github.com/rusandris/AttentionDynamics.jl.git
using AttentionDynamics
```

Or download a local copy, go to the directory and use the package in REPL:
```julia
] activate .
using AttentionDynamics
```

# References
Lorenz-Spreen, P., Mønsted, B.M., Hövel, P. et al. Accelerating dynamics of collective attention. Nat Commun 10, 1759 (2019). https://doi.org/10.1038/s41467-019-09311-w
