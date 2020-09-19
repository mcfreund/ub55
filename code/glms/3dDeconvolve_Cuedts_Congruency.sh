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
		dir_stimts=/stimts/${subject}/INPUT_DATA/Cuedts/${sessionsarray[$session_i]}
		dir_out=/out/${subject}/RESULTS/Cuedts/${sessionsarray[$session_i]}_CongruencyIncentive_EVENTS_censored_${runs[$run_i]}
		name_img=/img/${subject}/INPUT_DATA/Cuedts/${sessionsarray[$session_i]}/lpi_scale_blur4_tfMRI_Cuedts${sess}${runs[$run_i]}_${encoding_dir[$run_i]}.nii.gz

		## make result dir
		mkdir -p ${dir_out}

		## build xmat
		3dDeconvolve \
		-local_times \
		-x1D_stop \
		-GOFORIT 5 \
		-input ${name_img} \
		-polort A \
		-float \
		-censor ${dir_stimts}/movregs_FD_mask_run${runs[$run_i]}.txt \
		-num_stimts 6 \
		-stim_times_AM1 1 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_block_run${runs[$run_i]}.txt 'dmBLOCK(1)' -stim_label 1 block \
		-stim_times 2 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_blockONandOFF_run${runs[$run_i]}.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
		-stim_times 3 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_ConInc_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 3 ConInc \
		-stim_times 4 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_ConNoInc_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 4 ConNoInc \
		-stim_times 5 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_InConInc_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 5 InConInc \
		-stim_times 6 ${dir_stimts}/${subject}_Cuedts_${sessionsarray[$session_i]}_InConNoInc_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 6 InConNoInc \
		-ortvec ${dir_stimts}/motion_demean_${sessionsarray[$session_i]}_run${runs[$run_i]}.1D movregs \
		-x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
		-xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
		-nobucket

	done

done
