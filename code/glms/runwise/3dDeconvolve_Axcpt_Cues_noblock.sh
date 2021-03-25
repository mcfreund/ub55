#!/usr/bin/env bash

wd=$(pwd)

for session_i in ${!sessions[@]}; do

	for run_i in ${!runs[@]}; do
	

		## define paths and names
		sess=${sessions[$session_i]:0:3}  ## get short name
		sess=${sess^}  ## Namecase
		dir_stimts=${stimts}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}
		dir_out=${out}${subject}/RESULTS/Axcpt/${sessions[$session_i]}_Cues_EVENTS_censored_shifted_noblock_${runs[$run_i]}
		name_img=${img}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}/lpi_scale_tfMRI_Axcpt${sess}${runs[$run_i]}_${encoding_dir[$run_i]}_L.func.gii

		## make result dir
		mkdir -p ${dir_out}
		cd ${dir_out}
	
		## build xmat
        	/usr/local/pkg/afni_18/3dDeconvolve \
		-local_times \
		-force_TR 1.2 \
		-x1D_stop \
		-input ${name_img} \
		-polort A \
		-float \
		-censor ${dir_stimts}/movregs_FD_mask_run${runs[$run_i]}.txt \
		-num_stimts 7 \
		-stim_times 1 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_blockONandOFF_shifted_run${runs[$run_i]}.txt 'TENTzero(0,16.8,15)' -stim_label 1 blockONandOFF \
		-stim_times 2 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_AX_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 2 AX \
		-stim_times 3 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_AY_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 3 AY \
		-stim_times 4 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_Ang_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 4 Ang \
		-stim_times 5 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_BX_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 5 BX \
		-stim_times 6 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_BY_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 6 BY \
		-stim_times 7 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_Bng_shifted_run${runs[$run_i]}.txt 'TENTzero(0,21.6,19)' -stim_label 7 Bng \
		-ortvec ${dir_stimts}/Movement_Regressors_Axcpt${sess[$session_i]}${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
		-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
		-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
		-nobucket

	done

done

cd ${wd}
