#!/bin/bash

rs=( 4.25 3.75 )

for r in ${rs[*]}
	do
		echo $r
		julia --project precompute_trajectories_alpha_curve.jl $r > /dev/null 2> /dev/null &
	done
