#!/bin/bash

rs=( 6.5 6.25 5.75 5.5 5.25 4.75 4.0 )

for r in ${rs[*]}
	do
		echo $r
		julia --project precompute_trajectories_alpha_curve.jl $r > /dev/null 2> /dev/null &
	done
