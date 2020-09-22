#!/usr/bin/env bash

subjs=132017

for subj in ${subjs[@]}; do


	echo ${subj}
	
	
	singularity run \
#	-B /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/:/img \
#	-B /data/nil-external/ccp/freund/ub55/out/glms/:/stimts:ro \
#	-B /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_PREPROCESSED/${subj}/:/data:ro \
#	-B /data/nil-external/ccp/freund/ub55/out/glms/:/out \
	-B /data/nil-external/ccp/freund/ub55/code/glms/firstpass/:/scripts \
	/data/nil-bluearc/ccp-hcp/afni_analysis_aux.simg \
	--wave wave1 \
	--subject ${subj} \
	--session baseline \
	--task Axcpt Cuedts Stroop Stern \
	--origin /data/derivatives/fmriprep/ \
	--destination /out/ \
	--events /stimts/ \
	--pipeline fmriprep \
	--aux_analysis /scripts/ \
	--ncpus 1 > ./${subj}_runtime.log
	

done
