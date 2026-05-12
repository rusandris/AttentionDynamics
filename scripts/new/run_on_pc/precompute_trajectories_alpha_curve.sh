#!/bin/bash

rs=( 7.0 8.0 9.0 )

for r in ${rs[*]}
	do
		echo $r
		julia --project precompute_trajectories_alpha_curve.jl $r > /dev/null 2> /dev/null &
	done
