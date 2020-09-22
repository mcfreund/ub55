#!/usr/bin/env bash


## get vars

script_names=("3dDeconvolve_Axcpt_Cues")
sessions=baseline
#filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjects.txt"
#mapfile -t subjects < $filename
subjects=132017
runs=(1 2)
encoding_dir=(AP PA)


## directories

stimts=/data/nil-external/ccp/freund/ub55/out/glms/
out=/data/nil-external/ccp/freund/ub55/out/glms/
img=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/
scripts=/data/nil-external/ccp/freund/ub55/code/glms/


## loop over subjs


for subject in ${subjects[@]}
do
	echo ${subject}
	
	for script in ${script_names}
	do
	
		source ${script}
		
	done

done




for subject in ${subjects[@]}; do echo ${subject}; for script in ${script_names}; do source ${script}.sh;	done; done
