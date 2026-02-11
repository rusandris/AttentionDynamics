module AttentionDynamics
	using DynamicalSystems
	using SparseArrays
	using StatsBase
	using DataStructures
	using Colors
	using Setfield
	using BifurcationKit
	const BK=BifurcationKit
	using DifferentialEquations
	using LinearAlgebra
	using DelimitedFiles

	include("map.jl")
	include("produce_od.jl")
	include("evolve_system.jl")
	include("import_system.jl")
	include("count_limitcycles.jl")
	include("lyapunov.jl")
	include("find_periodic_orbit.jl")
	include("jacobian.jl")
	include("scaling_laws.jl")


	export map_to_integer,map_timeseries_to_integer,map_timeseries_to_binarystate,map_to_integer,map_to_binarystate,bifurcation_diagram,filter_repetitions,flip_randomstates
	export orbit_diagram,save_od_data
	export create_ds,create_integrator
	export count_limitcycles,find_limitcycle,findall_limcycles,findall_limcycles_many_systems
	export create_pinteg,lyapunov_curve,lyapunov_curves,lyapunov_diagram,replace_infs!
	export find_periodic_orbit,periodic_orbit_continuation,jacobian
	export period_scaling,amplitude_scaling
	
end # module AttentionDynamics
