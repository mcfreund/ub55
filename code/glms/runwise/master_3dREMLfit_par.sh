#!/usr/bin/env bash


## get vars

glms=(Cues CongruencySwitch ListLength Congruency Cues CongruencySwitch ListLength Congruency)
tasks=(Axcpt Cuedts Stern Stroop Axcpt Cuedts Stern Stroop)
suffices=(_shifted _shifted _shifted _shifted _shifted_noblock _shifted_noblock _shifted_noblock _shifted_noblock)
#suffices=(_shifted_noblock _shifted_noblock _shifted_noblock _shifted_noblock)
sessions=baseline
#filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
filename="/data/nil-external/ccp/freund/ub55/code/glms/runwise/missing_subjs.txt"
mapfile -t subjects < $filename
runs=(1 2)
encoding_dir=(AP PA)
hemis=(L R)
#hemis=(R)

#glms=(Cues CongruencySwitch ListLength Congruency fix-item)
#tasks=(Axcpt Cuedts Stern Stroop Stroop)
#suffices=(_shifted _shifted _shifted _shifted "")
#subjects=132017


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

