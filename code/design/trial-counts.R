library(here)
library(data.table)

source(here("code", "read-behav.R"))

## axcpt ----

## cue*target ax by bx ay
## response determined
## by run
## probe ID??

## basic design:

with(axcpt, table(interaction(subj, run), trial.type))
unique(axcpt$probe)
unique(axcpt$cue)
setdiff(unique(axcpt$probe), unique(axcpt$cue))  ## cue and probe letters shared, except for X (unique probe)


## probe
axcpt.probe <- as.matrix(with(axcpt, table(interaction(subj, run), probe)))
colSums(axcpt.probe)
colSums(axcpt.probe == 0)  ## all subjs have at least 1 of each probe in each run

axcpt.cue <- as.matrix(with(axcpt, table(interaction(subj, run), cue)))
colSums(axcpt.cue)
colSums(axcpt.cue == 0)  ## good number of each cue letter (except for #s, which are sparse)


## cue*probe
axcpt.cueprobe <- as.matrix(with(axcpt, table(interaction(subj, run), interaction(cue, probe))))
colSums(axcpt.cueprobe)
colSums(axcpt.cueprobe == 0)  ## of course every combination (of 100) is not represented.

## visual  i


