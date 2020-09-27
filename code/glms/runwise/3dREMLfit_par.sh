#!/usr/bin/env bash

# https://stackoverflow.com/questions/17621798/linux-process-in-background-stopped-in-jobs
# for writing stdin/stdout to file, for background running


wd=$(pwd)

function remlfit {
	
	local session_i=$1
	local run_i=$2
	
	## define paths and names
	sess=${sessions[$session_i]:0:3}  ## get short name
	sess=${sess^}  ## Namecase
	dir_stimts=${stimts}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}
	dir_out=${out}${subject}/RESULTS/${task}/${sessions[$session_i]}_${glm}_EVENTS_censored${suffix}_${runs[$run_i]}
	name_img=${img}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}/lpi_scale_blur4_tfMRI_${task}${sess}${runs[$run_i]}_${encoding_dir[$run_i]}.nii.gz

	cd ${dir_out}	

	/usr/local/pkg/afni_18/3dREMLfit \
	-matrix ${dir_out}/X_${runs[$run_i]}.xmat.1D \
	-input ${name_img} \
	-Rvar ${dir_out}/stats_var_${subject}_${runs[$run_i]}_REML.nii.gz \
	-Rbuck ${dir_out}/STATS_${subject}_${runs[$run_i]}_REML.nii.gz \
	-rwherr ${dir_out}/wherr_${subject}_${runs[$run_i]}_REML.nii.gz \
	-rerrts ${dir_out}/errts_${subject}_${runs[$run_i]}_REML.nii.gz \
	-GOFORIT \
	-fout \
	-tout \
	-nobout \
	-noFDR \
	-verb

}


for session_i in ${!sessions[@]}; do

	for run_i in ${!runs[@]}; do
		
		logpath=${out}${subject}/RESULTS/${task}/${sessions[$session_i]}_${glm}_EVENTS_censored${suffix}_${runs[$run_i]}
		remlfit ${session_i} ${run_i} < /dev/null > ${logpath}/runtime.log 2>&1 &

	done

done


cd ${wd}

