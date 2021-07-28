#!/usr/bin/env bash


tasks=("Axcpt" "Cuedts" "Stern" "Stroop")
name_out="aggressive1-catwherr"
name_invcov="baseline_aggressive1_EVENTS_censored_shifted_est-concat"

#${name_task_i}"-"${name_out}

for task_i in ${!tasks[@]}; do
	
	name_task_i=${tasks[$task_i]}
	
	./1_est_simil_rsa.R --name_invcov "invcov_"${name_task_i}"_"${name_invcov} --name_out ${name_task_i}"-"$name_out -p -n
	
done
