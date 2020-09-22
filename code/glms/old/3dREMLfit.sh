#!/usr/bin/env bash


runs=(1 2)
encoding_dir=(AP PA)
subject=${subjects}
sessionsarray=($sessions)


for session_i in ${!sessionsarray[@]}; do

	for run_i in ${!runs[@]}; do
	
		## define paths and names
		sess=${sessionsarray[$session_i]:0:3}  ## get short name
		sess=${sess^}  ## Namecase
		dir_stimts=/stimts/${subject}/INPUT_DATA/${task}/${sessionsarray[$session_i]}
		dir_out=/out/${subject}/RESULTS/${task}/${sessionsarray[$session_i]}_${glmname}_EVENTS_censored_${runs[$run_i]}
		name_img=/img/${subject}/INPUT_DATA/${task}/${sessionsarray[$session_i]}/lpi_scale_blur4_tfMRI_${task}${sess}${runs[$run_i]}_${encoding_dir[$run_i]}.nii.gz

		3dREMLfit \
		-matrix ${dir_out}/X_${runs[$run_i]}.xmat.1D \
		-GOFORIT 5 \
		-input ${name_img} \
		-Rvar ${dir_out}/stats_var_${subject}_${runs[$run_i]}_REML.nii.gz \
		-Rbuck ${dir_out}/STATS_${subject}_${runs[$run_i]}_REML.nii.gz \
		-rwherr ${dir_out}/wherr_${subject}_${runs[$run_i]}_REML.nii.gz \
		-rerrts ${dir_out}/errts_${subject}_${runs[$run_i]}_REML.nii.gz \
		-fout \
		-tout \
		-nobout \
		-noFDR \
		-verb

	done

done
