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

blocks.axcpt  <- axcpt.l %>%  lapply(write.blocks, dir.analysis = dir.analysis, by.run = TRUE)
blocks.axcpt  <- cuedts.l %>% lapply(write.blocks, dir.analysis = dir.analysis, by.run = TRUE)
blocks.stern  <- stern.l %>%  lapply(write.blocks, dir.analysis = dir.analysis, by.run = TRUE)
blocks.stroop <- stroop.l %>% lapply(write.blocks, dir.analysis = dir.analysis, by.run = TRUE)




## MOVREGS ----

to.split <- expand.grid(
  subj = subjs,
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  session = "baseline"
  )

movregs <- split.movregs(to.split, dir.to = dir.analysis)
sum(movregs$has.unexpected.nrow)
sum(movregs$is.missing.dir)




## AXCPT ----

## create arguments dataframe:
## (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)

args.axcpt.events <- expand.grid(
  var.level    = c("Ang", "AX", "AY", "Bng", "BX", "BY"),
  name.var     = "trial.type",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = TRUE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

args.axcpt.events

results.axcpt.events <- lapply(axcpt.l, write.events, .args = args.axcpt.events)
results.axcpt.events <- bind_rows(results.axcpt.events, .id = "subj.session")




## CUEDTS ----

args.cuedts.events <- expand.grid(
  var.level    = c("ConRepeat", "ConSwitch", "InConRepeat", "InConSwitch", "trial1"),
  name.var     = "trial.type.switch",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = TRUE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

args.cuedts.events

results.cuedts.events <- lapply(cuedts.l, write.events, .args = args.cuedts.events)
results.cuedts.events <- bind_rows(results.cuedts.events, .id = "subj.session")




## STERN ----

args.stern.events <- expand.grid(
  var.level    = c("LL5NN", "LL5NP", "LL5RN", "not5NN", "not5NP", "not5RN"),
  name.var     = "load01.trial.type",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset.shifted",
  by.run       = TRUE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

args.stern.events

results.stern.events <- lapply(stern.l, write.events, .args = args.stern.events)
results.stern.events <- bind_rows(results.stern.events, .id = "subj.session")




## STROOP ----

###### events
args.stroop.events <- expand.grid(
  var.level    = c("biasInCon", "biasCon", "PC50InCon", "PC50Con"),
  name.var     = "pc.trial.type",
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset.shifted",
  by.run       = TRUE,
  fname.suffix = "_shifted",
  stringsAsFactors = FALSE
)

args.stroop.events

results.stroop.events <- lapply(stroop.l, write.events, .args = args.stroop.events)
results.stroop.events <- bind_rows(results.stroop.events, .id = "subj.session")




## write records ----

## combine all results objects into single data.frame

n <- ls(pattern = "results")
l <- setNames(lapply(n, get), n)
l <- do.call(rbind, l)

fwrite(l, here("out", "glms", paste0("summary_write-events.csv")))

