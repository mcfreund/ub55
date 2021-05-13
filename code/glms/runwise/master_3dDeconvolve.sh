#!/usr/bin/env bash

### TODO:
# remove hemisphere loop---same glm, so unneccessary.
# make 3dDeconvolve parallel--- embed in function, call with forking.

## get vars

#glm_names=(Axcpt_Cues Cuedts_CongruencySwitch Stern_ListLength Stroop_Congruency Stroop_fix-item)
#glm_names=(Axcpt_Cues_noblock Cuedts_CongruencySwitch_noblock Stern_ListLength_noblock Stroop_Congruency_noblock)
glm_names=Cuedts_cueletnum
tasks=(Cuedts)
sessions=baseline
filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
mapfile -t subjects < $filename
#subjects=132017
runs=(1 2)
encoding_dir=(AP PA)

## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/ub55/code/glms/


## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm in ${glm_names[@]}; do
	
		source 3dDeconvolve_${glm}.sh

	done

done


