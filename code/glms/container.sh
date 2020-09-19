#!/usr/bin/env bash


filename="/data/nil-external/ccp/freund/ub55/in/ub55_subjs.txt"
mapfile -t subjs < $filename

for subj in ${subjs[@]}; do


	echo ${subj}
	
	
	singularity run \
	-B /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/:/img:ro \
	-B /data/nil-external/ccp/freund/ub55/glms/stimts/:/stimts:ro \
	-B /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_PREPROCESSED/${subj}/:/data:ro \
	-B /data/nil-external/ccp/freund/ub55/glms/:/out \
	-B /data/nil-external/ccp/freund/ub55/code/glms/:/scripts \
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
	--ncpus 4 > ./${subj}_runtime.log
	

done