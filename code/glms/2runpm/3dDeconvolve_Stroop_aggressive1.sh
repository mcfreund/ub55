#!/usr/bin/env bash

wd=$(pwd)

task=Stroop

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
        -stim_times 3 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blackBLACK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 3 blackBLACK \
        -stim_times 4 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blackGREEN_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 4 blackGREEN \
        -stim_times 5 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blackPINK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 5 blackPINK \
        -stim_times 6 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blackYELLOW_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 6 blackYELLOW \
        -stim_times 7 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blueBLUE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 7 blueBLUE \
        -stim_times 8 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_bluePURPLE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 8 bluePURPLE \
        -stim_times 9 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blueRED_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 9 blueRED \
        -stim_times 10 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_blueWHITE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 10 blueWHITE \
        -stim_times 11 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_greenBLACK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 11 greenBLACK \
        -stim_times 12 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_greenGREEN_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 12 greenGREEN \
        -stim_times 13 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_greenPINK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 13 greenPINK \
        -stim_times 14 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_greenYELLOW_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 14 greenYELLOW \
        -stim_times 15 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_pinkBLACK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 15 pinkBLACK \
        -stim_times 16 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_pinkGREEN_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 16 pinkGREEN \
        -stim_times 17 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_pinkPINK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 17 pinkPINK \
        -stim_times 18 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_pinkYELLOW_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 18 pinkYELLOW \
        -stim_times 19 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_purpleBLUE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 19 purpleBLUE \
        -stim_times 20 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_purplePURPLE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 20 purplePURPLE \
        -stim_times 21 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_purpleRED_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 21 purpleRED \
        -stim_times 22 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_purpleWHITE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 22 purpleWHITE \
        -stim_times 23 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_redBLUE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 23 redBLUE \
        -stim_times 24 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_redPURPLE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 24 redPURPLE \
        -stim_times 25 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_redRED_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 25 redRED \
        -stim_times 26 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_redWHITE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 26 redWHITE \
        -stim_times 27 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_whiteBLUE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 27 whiteBLUE \
        -stim_times 28 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_whitePURPLE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 28 whitePURPLE \
        -stim_times 29 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_whiteRED_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 29 whiteRED \
        -stim_times 30 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_whiteWHITE_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 30 whiteWHITE \
        -stim_times 31 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_yellowBLACK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 31 yellowBLACK \
        -stim_times 32 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_yellowGREEN_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 32 yellowGREEN \
        -stim_times 33 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_yellowPINK_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 33 yellowPINK \
        -stim_times 34 ${dir_stimts}/${subject}_${task}_${sessions[$session_i]}_yellowYELLOW_shifted.txt 'TENTzero(0,16.8,15)' -stim_label 34 yellowYELLOW \
        -ortvec ${dir_movreg}/motion_demean_baseline.1D movregs \
        -x1D ${dir_out}/X.xmat.1D \
        -xjpeg ${dir_out}/X.jpg \
        -nobucket

done


cd ${wd}
