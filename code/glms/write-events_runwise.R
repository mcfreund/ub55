## about ----
##
## writes events for run-wise GLMs.
##
## mike freund, 2020-04-04

library(here)
library(dplyr)
library(magrittr)
library(data.table)
library(mikeutils)

## functions

source(here("freund", "sub-subj-glms", "src", "runwise", "write-events_funs.R"))

## data

d <- read_dmcc_behav()
for (name.i in names(d)) assign(name.i, value = d[[name.i]])
rm(name.i, d)


## format ----
## 
## this section adapts columns within behav-and-event sheets to match the format of DMCC naming conventions.

axcpt$task  <- "Axcpt"
cuedts$task <- "Cuedts"
stern$task  <- "Stern"
stroop$task <- "Stroop"

axcpt$session  %<>% vapply(switch, character(1), "bas" = "baseline", "pro" = "proactive", "rea" = "reactive")
cuedts$session %<>% vapply(switch, character(1), "bas" = "baseline", "pro" = "proactive", "rea" = "reactive")
stern$session  %<>% vapply(switch, character(1), "bas" = "baseline", "pro" = "proactive", "rea" = "reactive")
stroop$session %<>% vapply(switch, character(1), "bas" = "baseline", "pro" = "proactive", "rea" = "reactive")

cuedts$trial.type %<>% vapply(switch, character(1), "c" = "Con", "i" = "InCon")
cuedts$incentive  %<>% vapply(switch, character(1), "nonincentive" = "NoInc", "incentive" = "Inc")
cuedts$incentive[cuedts$session == "baseline" & cuedts$target.color.orig == "green"] <- "Inc"
cuedts$switch[cuedts$switch == "stay"]   <- "Repeat"
cuedts$switch[cuedts$switch == "switch"] <- "Switch"
# substr(cuedts$switch, 1, 1) <- toupper(substr(cuedts$switch, 1, 1))

stroop$trial.type %<>% vapply(switch, character(1), "c" = "Con", "i" = "InCon")
is.pc50 <- stroop$pc == "pc50"
is.bias <- stroop$pc == "mi" | (stroop$pc == "mc" & stroop$session == "baseline")
is.buff <- stroop$pc == "mc" & stroop$session == "reactive"
stroop$pc[is.pc50]   <- "PC50"
stroop$pc[is.bias]   <- "bias"
stroop$pc[is.buff]   <- "buff"


## subjs to remove:
subj.incomplete  <- unique(subjsum$subj[subjsum$num.gii == 0 | subjsum$missing.rows | is.na(subjsum$missing.rows)])
subj.mb8         <- unique(subjsum$subj[subjsum$mb == "eight"])
subjsum.complete <- filter(subjsum, !subj %in% union(subj.incomplete, subj.mb8))
subj.complete    <- unique(subjsum.complete$subj)
subj.complete %<>% .[!. %in% c("155938", "300618")]  ## missing movregs

## select subjects for piloting:
# par(mfrow = c(1, 2))
# plot(subjsum.complete$mean.fd, subjsum.complete$per.censored)
# abline(h = 10, lty = "dashed", col = "firebrick")
# abline(v = 0.5, lty = "dashed", col = "firebrick")
# plot(subjsum$per.error)
# abline(h = 25, lty = "dashed", col = "firebrick")
# no.bueno <- which(subjsum$mean.fd > 0.5 | subjsum$per.censored > 10 | subjsum$per.error > 25)
# subj.no.bueno <- unique(subjsum$subj[no.bueno])
# subj.bueno <- subj.complete[!subj.complete %in% subj.no.bueno]
# subj.samp <- sample(subj.bueno, 5)
# subjsum[subjsum$subj %in% subj.samp, c("subj", "twin.pair")] %>% table


## prep dfs for loops
axcpt  %<>% filter(subj %in% subj.complete) %>% arrange(subj, session, trial.type)
cuedts %<>% filter(subj %in% subj.complete) %>% arrange(subj, session, trial.type)
stern  %<>% filter(subj %in% subj.complete) %>% arrange(subj, session, trial.type)
stroop %<>% filter(subj %in% subj.complete) %>% arrange(subj, session, trial.type)

dir.analysis <- "/data/nil-external/ccp/freund/sub-subj-glms/runwise_new"


## BLOCKS ----

axcpt.l <- split(axcpt, interaction(axcpt$subj, axcpt$session))
blocks.axcpt <- lapply(axcpt.l, write.blocks, dir.analysis = dir.analysis, by.run = TRUE)

cuedts.l <- split(cuedts, interaction(cuedts$subj, cuedts$session))
blocks.axcpt <- lapply(cuedts.l, write.blocks, dir.analysis = dir.analysis, by.run = TRUE)

stern.l <- split(stern, interaction(stern$subj, stern$session))
blocks.stern <- lapply(stern.l, write.blocks, dir.analysis = dir.analysis, by.run = TRUE)

stroop.l <- split(stroop, interaction(stroop$subj, stroop$session))
blocks.stroop <- lapply(stroop.l, write.blocks, dir.analysis = dir.analysis, by.run = TRUE)


## MOVREGS ----


to.split <- expand.grid(
  subj = subj.complete,
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  session = c("baseline", "proactive", "reactive")
  )

movregs <- split.movregs(to.split, dir.to = dir.analysis)
sum(movregs$has.unexpected.nrow)
sum(movregs$is.missing.dir)


## AXCPT ----

#### CREATE ARGUMENTS DATAFRAME --
####
#### (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)

###### buttons

args.axcpt.buttons <- expand.grid(
  var.level    = c("button1", "button2"),
  name.var     = "button",
  dir.analysis = dir.analysis,
  name.onset   = "rt",
  by.run       = TRUE,
  fname.suffix = "",
  stringsAsFactors = FALSE
)

args.axcpt.buttons

###### errors

args.axcpt.errors <- expand.grid(
  var.level    = "error",
  name.var     = "error",
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  fname.suffix = "",
  stringsAsFactors = FALSE
)

args.axcpt.errors

###### events

args.axcpt.events <- expand.grid(
  var.level    = c("Ang", "AX", "AY", "Bng", "BX", "BY"),
  name.var     = c("trial.type", "trial.type.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.axcpt.events$fname.suffix <- ifelse(grepl("correct", args.axcpt.events$name.var), "_correct", "")

args.axcpt.events


#### CREATE GROUPING VARS IN EVENTS DATAFRAME --
#### to identify which events to include

###### buttons
###### NB. buttons is different here (axcpt) because there are two button-presses for each trial.
###### thus, a separate events data.frame is needed for axcpt.buttons
axcpt.buttons <- axcpt %>%
  select(subj, task, session, run, trial.num, target.resp, cue.resp, time.target.rt, time.cue.rt) %>%
  as.data.table %>%
  melt(
    idvars = c("subj", "task", "session", "run", "trial.num"), 
    variable.name = "stim",
    measure.vars = list(rt = c("time.target.rt", "time.cue.rt"), resp = c("target.resp", "cue.resp"))
  ) %>% 
  mutate(
    button = ifelse(is.na(resp), NA, paste0("button", resp)),
    rt     = ifelse(is.na(resp), NA, rt)
    )

###### errors
axcpt$error <- ifelse(axcpt$trialacc == 0, "error", NA)

###### events
axcpt$trial.type.correct <- ifelse(axcpt$trialacc == 0, NA, axcpt$trial.type)



#### WRITE --

###### buttons
axcpt.buttons.l <- split(axcpt.buttons, list(axcpt.buttons$subj, axcpt.buttons$session))
results.axcpt.buttons <- lapply(axcpt.buttons.l, write.events, .args = args.axcpt.buttons)
results.axcpt.buttons <- bind_rows(results.axcpt.buttons, .id = "subj.session")

###### errors
axcpt.l <- split(axcpt, list(axcpt$subj, axcpt$session))
results.axcpt.errors <- lapply(axcpt.l, write.events, .args = args.axcpt.errors)
results.axcpt.errors <- bind_rows(results.axcpt.errors, .id = "subj.session")

###### events
results.axcpt.events <- lapply(axcpt.l, write.events, .args = args.axcpt.events)
results.axcpt.events <- bind_rows(results.axcpt.events, .id = "subj.session")



## CUEDTS ----

#### CREATE ARGUMENTS DATAFRAME --
####
#### (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)

###### buttons
args.cuedts.buttons <- args.axcpt.buttons
args.cuedts.buttons$name.onset <- "time.rt"

###### errors
args.cuedts.errors <- args.axcpt.errors
args.cuedts.errors

###### events
args.cuedts.trial.type.incentive <- expand.grid(
  var.level    = c("ConInc", "ConNoInc", "InConInc", "InConNoInc"),
  name.var     = c("trial.type.incentive", "trial.type.incentive.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.cuedts.switch.incentive <- expand.grid(
  var.level    = c("RepeatInc", "RepeatNoInc", "SwitchInc", "SwitchNoInc"),
  name.var     = c("switch.incentive", "switch.incentive.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.cuedts.trial.type.switch.incentive <- expand.grid(
  var.level    = c(
    "ConRepeatInc", "ConRepeatNoInc", "ConSwitchInc", "ConSwitchNoInc",
    "InConRepeatInc", "InConRepeatNoInc", "InConSwitchInc", "InConSwitchNoInc"
  ),
  name.var     = c("trial.type.switch.incentive", "trial.type.switch.incentive.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.cuedts.trial1 <- expand.grid(
  var.level    = "trial1",
  name.var     = c("trial1", "trial1.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)

args.cuedts.events <- rbind(
  args.cuedts.trial.type.incentive,
  args.cuedts.switch.incentive,
  args.cuedts.trial.type.switch.incentive,
  args.cuedts.trial1
)
args.cuedts.events$fname.suffix <- ifelse(grepl("correct", args.cuedts.events$name.var), "_correct", "")
args.cuedts.events


#### CREATE GROUPING VARS IN EVENTS DATAFRAME --
#### to identify which events to include

###### buttons
cuedts$button <- ifelse(is.na(cuedts$resp), NA, paste0("button", cuedts$resp))

###### errors
cuedts$error <- ifelse(cuedts$acc == 0, "error", NA)

###### events
cuedts$trial.type.incentive        <- paste0(cuedts$trial.type, cuedts$incentive)
cuedts$switch.incentive            <- ifelse(is.na(cuedts$switch), NA, paste0(cuedts$switch, cuedts$incentive))
cuedts$trial.type.switch.incentive <- ifelse(is.na(cuedts$switch), NA, paste0(cuedts$trial.type, cuedts$switch, cuedts$incentive))
# table(cuedts$trial.num, is.na(cuedts$switch))
cuedts$trial1                      <- ifelse(cuedts$trial.num %in% c(1, 19, 37), "trial1", NA)  ## beginning of blocks



cuedts$trial.type.incentive.correct        <- ifelse(cuedts$acc == 0, NA, cuedts$trial.type.incentive) 
cuedts$switch.incentive.correct            <- ifelse(cuedts$acc == 0, NA, cuedts$switch.incentive)
cuedts$trial.type.switch.incentive.correct <- ifelse(cuedts$acc == 0, NA, cuedts$trial.type.switch.incentive)
cuedts$trial1.correct                      <- ifelse(cuedts$acc == 0, NA, cuedts$trial1)

#### WRITE --

cuedts.l <- split(cuedts, list(cuedts$subj, cuedts$session))

###### buttons
results.cuedts.buttons <- lapply(cuedts.l, write.events, .args = args.cuedts.buttons)
results.cuedts.buttons <- bind_rows(results.cuedts.buttons, .id = "subj.session")


###### errors
results.cuedts.errors <- lapply(cuedts.l, write.events, .args = args.cuedts.errors)
results.cuedts.errors <- bind_rows(results.cuedts.errors, .id = "subj.session")

###### events
results.cuedts.events <- lapply(cuedts.l, write.events, .args = args.cuedts.events)
results.cuedts.events <- bind_rows(results.cuedts.events, .id = "subj.session")



## STERN ----

#### CREATE ARGUMENTS DATAFRAME --
####
#### (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)

###### buttons
args.stern.buttons <- args.cuedts.buttons
args.stern.buttons

###### errors
args.stern.errors <- args.axcpt.errors
args.stern.errors

###### events
args.stern.events <- expand.grid(
  var.level    = c("LL5NN", "LL5NP", "LL5RN", "not5NN", "not5NP", "not5RN"),
  name.var     = c("listlevel.trial.type", "listlevel.trial.type.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.cue.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.stern.events$fname.suffix <- ifelse(grepl("correct", args.stern.events$name.var), "_correct", "")
args.stern.events


#### CREATE GROUPING VARS IN EVENTS DATAFRAME --
#### to identify which events to include

###### buttons
stern$button <- ifelse(is.na(stern$resp), NA, paste0("button", stern$resp))

###### errors
stern$error <- ifelse(stern$acc == 0, "error", NA)

###### events
stern$listlevel.trial.type <- paste0(ifelse(stern$load == 5, "LL5", "not5"), stern$trial.type)
stern$listlevel.trial.type.correct <- ifelse(stern$acc == 0, NA, stern$listlevel.trial.type)


#### WRITE --

stern.l <- split(stern, list(stern$subj, stern$session))

###### buttons
results.stern.buttons <- lapply(stern.l, write.events, .args = args.stern.buttons)
results.stern.buttons <- bind_rows(results.stern.buttons, .id = "subj.session")

###### errors
results.stern.errors <- lapply(stern.l, write.events, .args = args.stern.errors)
results.stern.errors <- bind_rows(results.stern.errors, .id = "subj.session")

###### events
results.stern.events <- lapply(stern.l, write.events, .args = args.stern.events)
results.stern.events <- bind_rows(results.stern.events, .id = "subj.session")



## STROOP ----

#### CREATE ARGUMENTS DATAFRAME --
####
#### (args to create: var.level, name.var, fname.suffix, dir.analysis, name.onset, by.run)

###### errors
args.stroop.errors <- args.axcpt.errors
args.stroop.errors$name.onset <- "time.target.onset"

###### events
args.stroop.events <- expand.grid(
  var.level    = c("biasInCon", "biasCon", "PC50InCon", "PC50Con"),
  name.var     = c("pc.trial.type", "pc.trial.type.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.stroop.events$fname.suffix <- ifelse(grepl("correct", args.stroop.events$name.var), "_correct", "")
args.stroop.events

args.stroop.events.rea <- expand.grid(
  var.level    = "buffCon",
  name.var     = c("pc.trial.type", "pc.trial.type.correct"),
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset",
  by.run       = TRUE,
  stringsAsFactors = FALSE
)
args.stroop.events.rea$fname.suffix <- ifelse(grepl("correct", args.stroop.events.rea$name.var), "_correct", "")
args.stroop.events.rea

#### CREATE GROUPING VARS IN EVENTS DATAFRAME --
#### to identify which events to include

###### errors
stroop$error <- ifelse(stroop$acc == 0, "error", NA)

###### events
stroop$pc.trial.type <- paste0(stroop$pc, stroop$trial.type)
stroop$pc.trial.type.correct <- ifelse(stroop$acc == 0, NA, stroop$pc.trial.type)


#### WRITE --

stroop.l <- split(stroop, list(stroop$subj, stroop$session))

###### errors
results.stroop.errors <- lapply(stroop.l, write.events, .args = args.stroop.errors)
results.stroop.errors <- bind_rows(results.stroop.errors, .id = "subj.session")

###### events
results.stroop.events <- lapply(stroop.l, write.events, .args = args.stroop.events)
results.stroop.events <- bind_rows(results.stroop.events, .id = "subj.session")

stroop.l.rea <- stroop.l[-grep("baseline|proactive", names(stroop.l))]

results.stroop.events.rea <- lapply(stroop.l.rea, write.events, .args = args.stroop.events.rea)
results.stroop.events.rea <- bind_rows(results.stroop.events.rea, .id = "subj.session")

## write records ----

## combine all results objects into single data.frame

n <- ls(pattern = "results")
l <- setNames(lapply(n, get), n)
l <- do.call(rbind, l)

fwrite(l, here("freund", "sub-subj-glms", "out", "runwise", paste0("summary_write-events.csv")))

