dont source me!

library(dplyr)
library(data.table)
library(here)


ub55 <- c(
  "102008", "107321", "115825", "123117", "130114", "130518", "132017", "141422", "150423", "158136", "161832",
  "165032", "173738", "176845", "178950", "182436", "203418", "204319", "233326", "300618", "346945", "393550",
  "432332", "448347", "580650", "594156", "601127", "672756", "765864", "849971", "877168",
  "DMCC0472647", "DMCC1328342", "DMCC1624043", "DMCC3963378", "DMCC4854984", "DMCC5065053", "DMCC5195268", 
  "DMCC5244053", "DMCC5775387", "DMCC5820265", "DMCC6418065", "DMCC6627478", "DMCC6671683", "DMCC6705371", 
  "DMCC7297690", "DMCC7921988", "DMCC8050964", "DMCC8078683", "DMCC8214059", "DMCC8260571", "DMCC9441378",
  "DMCC9478705", "DMCC9850294", "DMCC9953810"
          )

dir.sheets <- "C:/Users/mcf/Box/DMCC_Phase2(HCP)/Preprocessed_Data/_wrangled"

stroop <- fread(file.path(dir.sheets, "dmcc2_behavior-and-events_stroop.csv"))[session == "bas" & subj %in% ub55]
axcpt  <- fread(file.path(dir.sheets, "dmcc2_behavior-and-events_axcpt.csv"))[session == "bas" & subj %in% ub55]
cuedts <- fread(file.path(dir.sheets, "dmcc2_behavior-and-events_cuedts.csv"))[session == "bas" & subj %in% ub55]
stern  <- fread(file.path(dir.sheets, "dmcc2_behavior-and-events_sternberg.csv"))[session == "bas" & subj %in% ub55]

setdiff(ub55, unique(stroop$subj))
setdiff(ub55, unique(axcpt$subj))
setdiff(ub55, unique(cuedts$subj))
setdiff(ub55, unique(stern$subj))

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
cuedts$switch[cuedts$switch == ""] <- NA
cuedts$trial.type.switch <- 
  ifelse(cuedts$trial.num %in% c(1, 19, 37), "trial1", paste0(cuedts$trial.type, cuedts$switch))

stern$load01 <- ifelse(stern$load == 5, "LL5", "not5")
stern$load01.trial.type <- paste0(stern$load01, stern$trial.type)


stroop$trial.type %<>% vapply(switch, character(1), "c" = "Con", "i" = "InCon")
is.pc50 <- stroop$pc == "pc50"
is.bias <- stroop$pc == "mi" | (stroop$pc == "mc" & stroop$session == "baseline")
is.buff <- stroop$pc == "mc" & stroop$session == "reactive"
stroop$pc[is.pc50]   <- "PC50"
stroop$pc[is.bias]   <- "bias"
stroop$pc[is.buff]   <- "buff"
stroop$pc.trial.type <- paste0(stroop$pc, stroop$trial.type)


fwrite(stroop, here("in", "ub55_stroop_behav-events.csv"))
fwrite(axcpt, here("in", "ub55_axcpt__behav-events.csv"))
fwrite(cuedts, here("in", "ub55_cuedts_behav-events.csv"))
fwrite(stern, here("in", "ub55_stern_behav-events.csv"))
fwrite(as.data.frame(ub55), here("in", "ub55_subjects.txt"), quote = FALSE, col.names = FALSE)
