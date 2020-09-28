#!/usr/bin/env bash


## get vars

glms=(fix-LSA)
tasks=(Stroop)
suffices=("")
sessions=baseline
#subjects=132017
filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
mapfile -t subjects < $filename
runs=(1 2)
encoding_dir=(AP PA)
hemis=(L R)


## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/ub55/code/

## fit


for subject in ${subjects[@]}; do

	echo ${subject}
	
	for glm_i in ${!glms[@]}; do
	
		for hemi in ${hemis[@]}; do

			glm=${glms[$glm_i]}
			task=${tasks[$glm_i]}
			suffix=${suffices[$glm_i]}
			
			source ${scripts}glms/runwise/3dREMLfit_par.sh
		
		done

	done

	wait

done

