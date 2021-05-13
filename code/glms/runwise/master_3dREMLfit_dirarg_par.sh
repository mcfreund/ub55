#!/usr/bin/env bash


## get vars

glms=(cueletnum)
tasks=(Cuedts)
suffices=("_shifted")
sessions=baseline
filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
#filename="/data/nil-external/ccp/freund/ub55/code/glms/runwise/missing_subjs.txt"
mapfile -t subjects < $filename
runs=(1 2)
encoding_dir=(AP PA)
hemis=(L R)
dirname="_EVENTS_censored"  ## e.g., _EVENTS_censored; in-between glm and suffix name

## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/


## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm_i in ${!glms[@]}; do
	
		glm=${glms[$glm_i]}
		task=${tasks[$glm_i]}
		suffix=${suffices[$glm_i]}
		
		echo ${task}_${glm}${suffix}

		source 3dREMLfit_dirarg_par.sh

	done

	wait

done

