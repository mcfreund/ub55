#!/usr/bin/env bash


## get vars

glms=(aggressive1 aggressive1 aggressive1 aggressive1)
tasks=(Axcpt Cuedts Stern Stroop)
suffices=(_shifted _shifted _shifted _shifted)
sessions=baseline
filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
mapfile -t subjects < $filename
hemis=(L R)


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

		source 3dREMLfit_par.sh

	done

	wait

done
