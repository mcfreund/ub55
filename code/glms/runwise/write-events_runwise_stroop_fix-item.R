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

stroop %<>% arrange(subj, session, trial.type)

## split

stroop.l <- stroop %>% split(interaction(.$subj, .$session))


## events

args.stroop.events <- expand.grid(
  var.level    = unique(stroop$item),
  name.var     = "item",
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset",
  by.run       = TRUE,
  fname.suffix = "",
  stringsAsFactors = FALSE
)

args.stroop.events

results.stroop.events <- lapply(stroop.l, write.events, .args = args.stroop.events)
results.stroop.events <- bind_rows(results.stroop.events, .id = "subj.session")

## write records

fwrite(results.stroop.events, here("out", "glms", paste0("summary_write-events_runwise_stroop_fix-item.csv")))
