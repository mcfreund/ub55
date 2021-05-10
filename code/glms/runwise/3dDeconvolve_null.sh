#!/usr/bin/env bash

wd=$(pwd)

for task_i in ${!tasks[@]}; do


	for session_i in ${!sessions[@]}; do

		for run_i in ${!runs[@]}; do
		

			## define paths and names
			sess=${sessions[$session_i]:0:3}  ## get short name
			sess=${sess^}  ## Namecase
			dir_stimts=${stimts}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}
			dir_out=${out}${subject}/RESULTS/${tasks[$task_i]}/${sessions[$session_i]}_null_${runs[$run_i]}
			name_img=${img}${subject}/INPUT_DATA/${tasks[$task_i]}/${sessions[$session_i]}/lpi_scale_tfMRI_${tasks[$task_i]}${sess}${runs[$run_i]}_${encoding_dir[$run_i]}_L.func.gii

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
			-num_stimts 0 \
			-ortvec ${dir_stimts}/Movement_Regressors_${tasks[$task_i]}${sess[$session_i]}${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
			-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
			-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
			-nobucket

		done

	done



done

cd ${wd}
