#!/usr/bin/env bash

wd=$(pwd)

for session_i in ${!sessions[@]}; do

	for run_i in ${!runs[@]}; do
	
		## define paths and names
		sess=${sessions[$session_i]:0:3}  ## get short name
		sess=${sess^}  ## Namecase
		dir_stimts=${stimts}${subject}/INPUT_DATA/Stroop/${sessions[$session_i]}
		dir_out=${out}${subject}/RESULTS/Stroop/${sessions[$session_i]}_fix-LSA_EVENTS_censored_${runs[$run_i]}
		name_img=${img}${subject}/INPUT_DATA/Stroop/${sessions[$session_i]}/lpi_scale_blur4_tfMRI_Stroop${sess}${runs[$run_i]}_${encoding_dir[$run_i]}.nii.gz

		## make result dir
		mkdir -p ${dir_out}
  		cd ${dir_out}


			## build xmat
            	/usr/local/pkg/afni_18/3dDeconvolve \
		-local_times \
		-x1D_stop \
		-input ${name_img} \
		-polort A \
		-float \
		-censor ${dir_stimts}/movregs_FD_mask_run${runs[$run_i]}.txt \
		-num_stimts 6 \
		-stim_times_AM1 1 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_block_run${runs[$run_i]}.txt 'dmBLOCK(1)' -stim_label 1 block \
		-stim_times_IM 2 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_blockONandOFF_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 2 blockONandOFF \
		-stim_times_IM 3 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_PC50Con_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 3 PC50Con \
		-stim_times_IM 4 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_PC50InCon_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 4 PC50InCon \
		-stim_times_IM 5 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_biasCon_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 5 biasCon \
		-stim_times_IM 6 ${dir_stimts}/${subject}_Stroop_${sessions[$session_i]}_biasInCon_run${runs[$run_i]}.txt 'BLOCK(1,1)' -stim_label 6 biasInCon \
            	-ortvec ${dir_stimts}/Movement_Regressors_Stroop${sess[$session_i]}${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
		-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
		-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
		-nobucket

	done

done

cd ${wd}
