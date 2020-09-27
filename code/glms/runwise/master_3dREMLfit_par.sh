#!/usr/bin/env bash


## get vars

glms=(Cues CongruencySwitch ListLength Congruency fix-item)
tasks=(Axcpt Cuedts Stern Stroop Stroop)
suffices=(_shifted _shifted _shifted _shifted "")
sessions=baseline
subjects=132017
runs=(1 2)
encoding_dir=(AP PA)


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
		
		source 3dREMLfit_par.sh

	done

done


