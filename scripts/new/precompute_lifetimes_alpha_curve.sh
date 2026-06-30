#!/bin/bash

rdirs=( B24-r3.75 B24-r4.0 B24-r4.25 B24-r5.25 B24-r5.5 B24-r5.75 B24-r6.25 B24-r6.5 )

for rdir in ${rdirs[*]}
        do
        	julia --project chart_stats_master_N300.jl $rdir > /dev/null &
	done  
