#!/usr/bin/env bash

subjs=132017

for subj in ${subjs[@]}; do


	echo ${subj}
	
	
	singularity run \
	-B /data/nil-external/ccp/freund/ub55/out/glms/:/stimts:ro \
	-B /data/nil-external/ccp/freund/ub55/out/glms/:/out \
	-B /data/nil-external/ccp/freund/ub55/code/glms/firstpass/:/scripts \
	/data/nil-bluearc/ccp-hcp/afni_analysis.simg \
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
