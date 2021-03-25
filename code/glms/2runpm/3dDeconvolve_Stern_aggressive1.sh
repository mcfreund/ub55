#!/usr/bin/env bash

wd=$(pwd)

task=Stern

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
        -num_stimts 13 \
        -stim_times_AM1 1 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_block.txt 'dmBLOCK(1)' -stim_label 1 block \
        -stim_times 2 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blockONandOFF_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
        -stim_times 3 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load5_NN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 3 load5_NN \
        -stim_times 4 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load6_NN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 4 load6_NN \
        -stim_times 5 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load7_NN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 5 load7_NN \
        -stim_times 6 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load8_NN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 6 load8_NN \
        -stim_times 7 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load5_NP_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 7 load5_NP \
        -stim_times 8 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load6_NP_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 8 load6_NP \
        -stim_times 9 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load7_NP_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 9 load7_NP \
        -stim_times 10 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load8_NP_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 10 load8_NP \
        -stim_times 11 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load5_RN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 11 load5_RN \
        -stim_times 12 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load7_RN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 12 load7_RN \
        -stim_times 13 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_load8_RN_shifted.txt 'TENTzero(0,26.4,23)' -stim_label 13 load8_RN \
        -ortvec ${dir_movreg}/motion_demean_baseline.1D movregs \
        -x1D ${dir_out}/X.xmat.1D \
        -xjpeg ${dir_out}/X.jpg \
        -nobucket

done


cd ${wd}
