#!/usr/bin/env bash

## get vars

glm_names=(fix-LSA)
tasks=(Stroop)
suffices=("")
sessions=baseline
#filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
#mapfile -t subjects < $filename
subjects=132017
runs=(1 2)
encoding_dir=(AP PA)
hemis=(L R)

## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/ub55/code/glms/


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm_i in ${!glm_names[@]}; do
	
		wd=$(pwd)

		for session_i in ${!sessions[@]}; do

			for run_i in ${!runs[@]}; do

				for hemi in ${hemis[@]}; do
				
					dir_out=${out}${subject}/RESULTS/${tasks[glm_i]}/${sessions[$session_i]}_${glm_names[$glm_i]}_EVENTS_censored${suffices[glm_i]}_${runs[$run_i]}/
					cd ${dir_out}
					atlas=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/Schaefer2018_400Parcels_7Networks_order_10K_${hemi}.label.gii
					stats=STATS_${subject}_${runs[$run_i]}_${hemi}_REML.func.gii

					3dROIstats -mask $atlas $stats > ${dir_out}roistats_schaefer400-07_${subject}_${runs[run_i]}_${hemi}.txt

				done  ## hemi

			done  ## run

			cd ${wd}

		done  ## session

	done  ## glm

done  ## subj


