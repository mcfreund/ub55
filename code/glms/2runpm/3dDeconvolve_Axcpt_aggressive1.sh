#!/usr/bin/env bash

wd=$(pwd)

for session_i in ${!sessions[@]}; do

	## define paths and names
	sess=${sessions[$session_i]:0:3}  ## get short name
	sess=${sess^}  ## Namecase
	dir_stimts=${stimts}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}
	dir_out=${out}${subject}/RESULTS/Axcpt/${sessions[$session_i]}_aggressive1_EVENTS_censored_shifted
	name_img1=${img}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}/lpi_scale_tfMRI_Axcpt${sess}1_AP_L.func.gii
	name_img2=${img}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}/lpi_scale_tfMRI_Axcpt${sess}2_PA_L.func.gii
	dir_movreg=${img}${subject}/INPUT_DATA/Axcpt/${sessions[$session_i]}

	## make result dir
	mkdir -p ${dir_out}
	cd ${dir_out}

	## build xmat
	/usr/local/pkg/afni_18/3dDeconvolve \
	-local_times \
	-force_TR 1.2 \
	-x1D_stop \
	-input ${name_img1} ${name_img2} \
	-polort A \
	-float \
	-censor ${dir_movreg}/movregs_FD_mask.txt \
	-num_stimts 32 \
	-stim_times_AM1 1 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_block.txt 'dmBLOCK(1)' -stim_label 1 block \
	-stim_times 2 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_blockONandOFF_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
	-stim_times 3 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_A_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 3 A_ng \
	-stim_times 4 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_C_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 4 C_ng \
	-stim_times 5 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_D_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 5 D_ng \
	-stim_times 6 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_F_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 6 F_ng \
	-stim_times 7 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_H_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 7 H_ng \
	-stim_times 8 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_M_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 8 M_ng \
	-stim_times 9 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_N_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 9 N_ng \
	-stim_times 10 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_P_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 10 P_ng \
	-stim_times 11 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_T_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 11 T_ng \
	-stim_times 12 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_U_ng_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 12 U_ng \
	-stim_times 13 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_A_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 13 A_X \
	-stim_times 14 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_C_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 14 C_X \
	-stim_times 15 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_D_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 15 D_X \
	-stim_times 16 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_F_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 16 F_X \
	-stim_times 17 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_H_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 17 H_X \
	-stim_times 18 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_M_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 18 M_X \
	-stim_times 19 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_N_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 19 N_X \
	-stim_times 20 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_P_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 20 P_X \
	-stim_times 21 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_T_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 21 T_X \
	-stim_times 22 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_U_X_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 22 U_X \
	-stim_times 23 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_A_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 23 A_Y \
	-stim_times 24 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_C_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 24 C_Y \
	-stim_times 25 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_D_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 25 D_Y \
	-stim_times 26 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_F_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 26 F_Y \
	-stim_times 27 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_H_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 27 H_Y \
	-stim_times 28 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_M_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 28 M_Y \
	-stim_times 29 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_N_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 29 N_Y \
	-stim_times 30 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_P_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 30 P_Y \
	-stim_times 31 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_T_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 31 T_Y \
	-stim_times 32 ${dir_stimts}/${subject}_Axcpt_${sessions[$session_i]}_U_Y_shifted.txt 'TENTzero(0,21.6,19)' -stim_label 32 U_Y \
	-ortvec ${dir_movreg}/motion_demean_baseline.1D movregs \
	-x1D ${dir_out}/X.xmat.1D \
	-xjpeg ${dir_out}/X.jpg \
	-nobucket

done

cd ${wd}
