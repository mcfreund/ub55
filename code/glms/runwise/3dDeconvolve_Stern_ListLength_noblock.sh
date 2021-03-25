#!/usr/bin/env bash

wd=$(pwd)

for session_i in ${!sessions[@]}; do

	for run_i in ${!runs[@]}; do
	
		## define paths and names
		sess=${sessions[$session_i]:0:3}  ## get short name
		sess=${sess^}  ## Namecase
		dir_stimts=${stimts}${subject}/INPUT_DATA/Stern/${sessions[$session_i]}
		dir_out=${out}${subject}/RESULTS/Stern/${sessions[$session_i]}_ListLength_EVENTS_censored_shifted_noblock_${runs[$run_i]}
		name_img=${img}${subject}/INPUT_DATA/Stern/${sessions[$session_i]}/lpi_scale_tfMRI_Stern${sess}${runs[$run_i]}_${encoding_dir[$run_i]}_L.func.gii

		## make result dir
		mkdir -p ${dir_out}
		cd ${dir_out}

		## build xmat
	        /usr/local/pkg/afni_18/3dDeconvolve \
		-local_times \
		-x1D_stop \
		-force_TR 1.2 \
		-input ${name_img} \
		-polort A \
		-float \
		-censor ${dir_stimts}/movregs_FD_mask_run${runs[$run_i]}.txt \
		-num_stimts 7 \
		-stim_times 1 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_blockONandOFF_shifted_run${runs[$run_i]}.txt 'TENTzero(0,16.8,15)' -stim_label 1 blockONandOFF \
		-stim_times 2 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_LL5NP_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 2 LL5NP \
		-stim_times 3 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_LL5NN_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 3 LL5NN \
		-stim_times 4 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_LL5RN_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 4 LL5RN \
		-stim_times 5 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_not5NP_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 5 not5NP \
		-stim_times 6 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_not5NN_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 6 not5NN \
		-stim_times 7 ${dir_stimts}/${subject}_Stern_${sessions[$session_i]}_not5RN_shifted_run${runs[$run_i]}.txt 'TENTzero(0,26.4,23)' -stim_label 7 not5RN \
	        -ortvec ${dir_stimts}/Movement_Regressors_Stern${sess[$session_i]}${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
		-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
		-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
		-nobucket
	
	done

done

cd ${wd}
