## about ----
##
## writes events for run-wise GLMs.
##
## mike freund, 2020-04-04
## adapted for ub55 2020-08-20

library(here)
library(dplyr)
library(magrittr)
library(data.table)
library(mikeutils)

source(here("code", "glms", "write-events_funs.R"))
source(here("code", "read-behav.R"))

dir.analysis <- here("out", "glms")



## format ----

## prep dfs for loops

axcpt  %<>% arrange(subj, session, trial.type)
cuedts %<>% arrange(subj, session, trial.type)
stern  %<>% arrange(subj, session, trial.type)
stroop %<>% arrange(subj, session, trial.type)

## shift events back in time to beginning of first TR of event

axcpt$time.cue.onset.shifted     <- axcpt$time.cue.onset - axcpt$time.cue.onset %% 1.2
cuedts$time.cue.onset.shifted    <- cuedts$time.cue.onset - cuedts$time.cue.onset %% 1.2
stern$time.cue.onset.shifted     <- stern$time.cue.onset - stern$time.cue.onset %% 1.2
stroop$time.target.onset.shifted <- stroop$time.target.onset - stroop$time.target.onset %% 1.2

## split

axcpt.l  <- axcpt %>%  split(interaction(.$subj, .$session))
cuedts.l <- cuedts %>% split(interaction(.$subj, .$session))
stern.l  <- stern %>%  split(interaction(.$subj, .$session))
stroop.l <- stroop %>% split(interaction(.$subj, .$session))




## BLOCKS ----

lapply(
  list(axcpt = axcpt, cuedts = cuedts, stern = stern, stroop = stroop),
  function(x) {
    
    block1on <- c(unique(x$time.block1.on - x$time.block1.on %% 1.2), range(unique(x$time.block1.on)))
    block2on <- c(unique(x$time.block2.on - x$time.block2.on %% 1.2), range(unique(x$time.block2.on)))
    block3on <- c(unique(x$time.block3.on - x$time.block3.on %% 1.2), range(unique(x$time.block3.on)))
    block1off <- c(unique(x$time.block1.off - x$time.block1.off %% 1.2), range(unique(x$time.block1.off)))
    block2off <- c(unique(x$time.block2.off - x$time.block2.off %% 1.2), range(unique(x$time.block2.off)))
    block3off <- c(unique(x$time.block3.off - x$time.block3.off %% 1.2), range(unique(x$time.block3.off)))
    
    b <- rbind(block1on, block2on, block3on, block1off, block2off, block3off)
    colnames(b) <- c("shift", "min", "max")
    
    b
    
  }
)

## TODO: 
## - block w/ single duration param per task, begin/end at actual block on/off times (get mean duration values)
## - 1trpk transient begin/end at shifted times
 
blocks.axcpt  <- axcpt.l %>%  lapply(write.blocks, dir.analysis = dir.analysis, by.run = FALSE)
blocks.axcpt  <- cuedts.l %>% lapply(write.blocks, dir.analysis = dir.analysis, by.run = FALSE)
blocks.stern  <- stern.l %>%  lapply(write.blocks, dir.analysis = dir.analysis, by.run = FALSE)
blocks.stroop <- stroop.l %>% lapply(write.blocks, dir.analysis = dir.analysis, by.run = FALSE)


blocks.shifted.axcpt  <- axcpt.l %>%  lapply(write.blocks.shifted, dir.analysis = dir.analysis, by.run = FALSE)
blocks.shifted.axcpt  <- cuedts.l %>% lapply(write.blocks.shifted, dir.analysis = dir.analysis, by.run = FALSE)
blocks.shifted.stern  <- stern.l %>%  lapply(write.blocks.shifted, dir.analysis = dir.analysis, by.run = FALSE)
blocks.shifted.stroop <- stroop.l %>% lapply(write.blocks.shifted, dir.analysis = dir.analysis, by.run = FALSE)


## MOVREGS ----
## will need to lookup movregs from stimtimes files on nil-bluearc



## AXCPT ----

## create arguments dataframe:
## (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)


for (ii in seq_along(axcpt.l)) {
  axcpt.l[[ii]]$cue_probetype <- paste0(axcpt.l[[ii]]$cue, "_", gsub("^.", "", axcpt.l[[ii]]$trial.type))
}

# unique(unlist(lapply(axcpt.l, function(x) unique(x$cue_probetype))))
cue_probetype <- combo_paste(c("A", "C", "D", "F", "H", "M", "N", "P", "T", "U"), c("ng", "X", "Y"))

args.axcpt.events <- expand.grid(
  var.level    = cue_probetype,  ## see trialcounts.html
  name.var     = "cue_probetype",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = FALSE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

## generate only shifted
# args.axcpt.events <- bind_rows(
#   args.axcpt.events,
#   args.axcpt.events %>% mutate(name.onset = gsub(".shifted", "", name.onset), fname.suffix = "")
# )

args.axcpt.events

results.axcpt.events <- lapply(axcpt.l, write.events, .args = args.axcpt.events)
results.axcpt.events <- bind_rows(results.axcpt.events, .id = "subj.session")

hist(as.numeric(results.axcpt.events$n.events))
plot(as.numeric(results.axcpt.events$n.events))


## CUEDTS ----

# lapply(cuedts.l, function(x) unique(paste0(x$cue, "_", x$stimuli))) %>% lapply(length)
cue_probeid <- combo_paste(c("l", "n"), unique(unlist(lapply(cuedts.l, function(x) unique(x$stimuli)))))

for (ii in seq_along(cuedts.l)) {
  cuedts.l[[ii]]$cue_probeid <- paste0(cuedts.l[[ii]]$cue, "_", cuedts.l[[ii]]$stimuli)
}


args.cuedts.events <- expand.grid(
  var.level    = cue_probeid,
  name.var     = "cue_probeid",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = FALSE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

# args.cuedts.events <- bind_rows(
#   args.cuedts.events,
#   args.cuedts.events %>% mutate(name.onset = gsub(".shifted", "", name.onset), fname.suffix = "")
# )

args.cuedts.events

results.cuedts.events <- lapply(cuedts.l, write.events, .args = args.cuedts.events)
results.cuedts.events <- bind_rows(results.cuedts.events, .id = "subj.session")
hist(as.numeric(results.cuedts.events$n.events))
plot(as.numeric(results.cuedts.events$n.events))



## STERN ----


for (ii in seq_along(stern.l)) {
  stern.l[[ii]]$load_trialtype <- paste0("load", stern.l[[ii]]$load, "_", stern.l[[ii]]$trial.type)
}


args.stern.events <- expand.grid(
  var.level    = load_trialtype,
  name.var     = "load_trialtype",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = FALSE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

# args.stern.events <- bind_rows(
#   args.stern.events,
#   args.stern.events %>% mutate(name.onset = gsub(".shifted", "", name.onset), fname.suffix = "")
# )

args.stern.events

results.stern.events <- lapply(stern.l, write.events, .args = args.stern.events)
results.stern.events <- bind_rows(results.stern.events, .id = "subj.session")
hist(as.numeric(results.stern.events$n.events))
plot(as.numeric(results.stern.events$n.events))



## STROOP ----

stroop_stimuli <- unique(unlist(lapply(stroop.l, function(x) unique(x$item))))

###### events

args.stroop.events <- expand.grid(
  var.level    = stroop_stimuli,
  name.var     = "item",
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset.shifted",
  by.run       = FALSE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

args.stroop.events <- bind_rows(
  args.stroop.events,
  args.stroop.events %>% mutate(name.onset = gsub(".shifted", "", name.onset), fname.suffix = "")
)

results.stroop.events <- lapply(stroop.l, write.events, .args = args.stroop.events)
results.stroop.events <- bind_rows(results.stroop.events, .id = "subj.session")
hist(as.numeric(results.stroop.events$n.events))
plot(as.numeric(results.stroop.events$n.events))



## write records ----

## combine all results objects into single data.frame

n <- ls(pattern = "results")
l <- setNames(lapply(n, get), n)
l <- do.call(rbind, l)

fwrite(l, here("out", "glms", paste0("summary_write-events_2runpm_agressive1.csv")))

