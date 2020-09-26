#!/usr/bin/env bash



## get vars

glm_names=(Stroop_fix-LSA)
sessions=baseline
filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
mapfile -t subjects < $filename
runs=(1 2)
encoding_dir=(AP PA)


## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/


## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm in ${glm_names[@]}; do
	
		source 3dDeconvolve_${glm}.sh

	done

done


