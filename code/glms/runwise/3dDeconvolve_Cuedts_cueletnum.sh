#!/usr/bin/env bash

wd=$(pwd)

for session_i in ${!sessions[@]}; do

    for run_i in ${!runs[@]}; do
	

        ## define paths and names
        sess=${sessions[$session_i]:0:3}  ## get short name
        sess=${sess^}  ## Namecase
        dir_stimts=${stimts}${subject}/INPUT_DATA/Cuedts/${sessions[$session_i]}
        dir_out=${out}${subject}/RESULTS/Cuedts/${sessions[$session_i]}_cueletnum_EVENTS_censored_shifted_${runs[$run_i]}
        name_img=${img}${subject}/INPUT_DATA/Cuedts/${sessions[$session_i]}/lpi_scale_tfMRI_Cuedts${sess}${runs[$run_i]}_${encoding_dir[$run_i]}_L.func.gii

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
        -num_stimts 10 \
        -stim_times_AM1 1 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_block_run${runs[$run_i]}.txt 'dmBLOCK(1)' -stim_label 1 block \
        -stim_times 2 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_blockONandOFF_shifted_run${runs[$run_i]}.txt 'TENTzero(0,16.8,15)' -stim_label 2 blockONandOFF \
        -stim_times 3 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_l_o_c_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 3 let_odd_con \
        -stim_times 4 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_l_o_v_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 4 let_odd_vow \
        -stim_times 5 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_l_e_c_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 5 let_eve_con \
        -stim_times 6 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_l_e_v_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 6 let_eve_vow \
        -stim_times 7 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_n_o_c_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 7 num_odd_con \
        -stim_times 8 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_n_o_v_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 8 num_odd_vow \
        -stim_times 9 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_n_e_c_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 9 num_eve_con \
        -stim_times 10 ${dir_stimts}/${subject}_Cuedts_${sessions[$session_i]}_n_e_v_shifted_run${runs[$run_i]}.txt 'TENTzero(0,24,21)' -stim_label 10 num_eve_vow \
        -ortvec ${dir_stimts}/Movement_Regressors_Cuedts${sess[$session_i]}${runs[$run_i]}_${encoding_dir[$run_i]}.1D movregs \
        -x1D ${dir_out}/X_${runs[$run_i]}.xmat.1D \
        -xjpeg ${dir_out}/X_${runs[$run_i]}.jpg \
        -nobucket
		
    done

done

cd ${wd}
