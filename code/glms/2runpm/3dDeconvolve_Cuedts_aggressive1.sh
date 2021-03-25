#!/usr/bin/env bash

wd=$(pwd)

task=Cuedts

for session_i in ${!sessions[@]}; do

        ## define paths and names
        sess=${sessions[$session_i]:0:3}  ## get short name
        sess=${sess^}  ## Namecase
        dir_stimts=${stimts}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}
        dir_out=${out}${subject}/RESULTS/${task}/${sessions[$session_i]}_aggressive1_EVENTS_censored_shifted
        name_img1=${img}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}/lpi_scale_tfMRI_${task}${sess}1_AP_L.func.gii
        name_img2=${img}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}/lpi_scale_tfMRI_${task}${sess}2_PA_L.func.gii
        dir_movreg=${img}${subject}/INPUT_DATA/${task}/${sessions[$session_i]}

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
        -num_stimts 34 \
        -stim_times_AM1 1 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_block.txt 'dmBLOCK(1)' -stim_label 1 block \
        -stim_times 2 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blockONandOFF_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
        -stim_times 3 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_6I_shifted.txt 'TENTzero(0,24,21)' -stim_label 3 l_6I \
        -stim_times 4 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_I6_shifted.txt 'TENTzero(0,24,21)' -stim_label 4 l_I6 \
        -stim_times 5 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_A2_shifted.txt 'TENTzero(0,24,21)' -stim_label 5 l_A2 \
        -stim_times 6 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_B1_shifted.txt 'TENTzero(0,24,21)' -stim_label 6 l_B1 \
        -stim_times 7 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_2A_shifted.txt 'TENTzero(0,24,21)' -stim_label 7 l_2A \
        -stim_times 8 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_1B_shifted.txt 'TENTzero(0,24,21)' -stim_label 8 l_1B \
        -stim_times 9 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_5H_shifted.txt 'TENTzero(0,24,21)' -stim_label 9 l_5H \
        -stim_times 10 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_H5_shifted.txt 'TENTzero(0,24,21)' -stim_label 10 l_H5 \
        -stim_times 11 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_B2_shifted.txt 'TENTzero(0,24,21)' -stim_label 11 l_B2 \
        -stim_times 12 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_D4_shifted.txt 'TENTzero(0,24,21)' -stim_label 12 l_D4 \
        -stim_times 13 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_4D_shifted.txt 'TENTzero(0,24,21)' -stim_label 13 l_4D \
        -stim_times 14 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_1A_shifted.txt 'TENTzero(0,24,21)' -stim_label 14 l_1A \
        -stim_times 15 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_A1_shifted.txt 'TENTzero(0,24,21)' -stim_label 15 l_A1 \
        -stim_times 16 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_E3_shifted.txt 'TENTzero(0,24,21)' -stim_label 16 l_E3 \
        -stim_times 17 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_3E_shifted.txt 'TENTzero(0,24,21)' -stim_label 17 l_3E \
        -stim_times 18 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_l_2B_shifted.txt 'TENTzero(0,24,21)' -stim_label 18 l_2B \
        -stim_times 19 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_6I_shifted.txt 'TENTzero(0,24,21)' -stim_label 19 n_6I \
        -stim_times 20 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_I6_shifted.txt 'TENTzero(0,24,21)' -stim_label 20 n_I6 \
        -stim_times 21 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_A2_shifted.txt 'TENTzero(0,24,21)' -stim_label 21 n_A2 \
        -stim_times 22 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_2A_shifted.txt 'TENTzero(0,24,21)' -stim_label 22 n_2A \
        -stim_times 23 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_B1_shifted.txt 'TENTzero(0,24,21)' -stim_label 23 n_B1 \
        -stim_times 24 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_1B_shifted.txt 'TENTzero(0,24,21)' -stim_label 24 n_1B \
        -stim_times 25 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_5H_shifted.txt 'TENTzero(0,24,21)' -stim_label 25 n_5H \
        -stim_times 26 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_H5_shifted.txt 'TENTzero(0,24,21)' -stim_label 26 n_H5 \
        -stim_times 27 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_B2_shifted.txt 'TENTzero(0,24,21)' -stim_label 27 n_B2 \
        -stim_times 28 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_D4_shifted.txt 'TENTzero(0,24,21)' -stim_label 28 n_D4 \
        -stim_times 29 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_4D_shifted.txt 'TENTzero(0,24,21)' -stim_label 29 n_4D \
        -stim_times 30 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_1A_shifted.txt 'TENTzero(0,24,21)' -stim_label 30 n_1A \
        -stim_times 31 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_A1_shifted.txt 'TENTzero(0,24,21)' -stim_label 31 n_A1 \
        -stim_times 32 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_E3_shifted.txt 'TENTzero(0,24,21)' -stim_label 32 n_E3 \
        -stim_times 33 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_2B_shifted.txt 'TENTzero(0,24,21)' -stim_label 33 n_2B \
        -stim_times 34 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_n_3E_shifted.txt 'TENTzero(0,24,21)' -stim_label 34 n_3E \
        -ortvec ${dir_movreg}/motion_demean_baseline.1D movregs \
        -x1D ${dir_out}/X.xmat.1D \
        -xjpeg ${dir_out}/X.jpg \
        -nobucket

done

cd ${wd}
