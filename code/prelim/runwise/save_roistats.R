

library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(magrittr)
library(gifti)
library(cifti)
library(abind)
library(mikeutils)
library(progress)

source(here("code", "glms", "write-events_funs.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))



## input: subjs, task, session, glmname
## output: RDS files

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop", "Stroop"),
  name.glm = c(
    "baseline_Cues_EVENTS_censored_shifted", 
    # "CongruencyIncentive_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted",
    "baseline_fix-item_EVENTS_censored"
  )
)
glminfo <- as.data.table(glminfo)
# dirs <- full_join(dirs, glminfo, by = "task")
# dirs$name.glm <- paste0(dirs$session, "_", dirs$name.glm)
# dirs$run1 <- file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, combo_paste(dirs$name.glm, c("1", "2"), c("L", "R")), "X_1.xmat.1D")
# 
# combopaste(
#   file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, dirs$name.glm), c("_1", "_2"), "/", 
#   "roistats_schaefer400-07_", "_1_L.txt"
#   ,
#   ".txt"
# )
# 
# 
# dirs <- pivot_longer(dirs, cols = c("run1", "run2"), names_to = "run", values_to = "fname.xmat")
# dirs$exists.xmat <- file.exists(dirs$fname.xmat)
# 
# dirs$fname.censor <- dirs$fname.xmat %>% gsub("RESULTS", "INPUT_DATA", .) %>% gsub("/baseline_.*", "/baseline/movregs_FD_mask_", .)
# dirs$fname.censor <- paste0(dirs$fname.censor, dirs$run, ".txt")
# 
# dirs <- filter(dirs, !subj %in% "432332")  ## for now, remove those that don't exist
# 
# dirs <- as.data.table(dirs)  ## for fast extracting


# file.path(dir.analysis, "132017", "RESULTS", "Axcpt", "baseline_Cues_EVENTS_censored_shifted_1", "roistats_schaefer400-07_132017_1_L.txt")

## loop over tasks, subjs, runs; write 3droistats to results

res <- vector("list", nrow(glminfo) * length(subjs) * 4)
names(res) <- combo_paste(paste0(glminfo$task, glminfo$name.glm), subjs, 1:2, c("L", "R"))

for (task.i in seq_along(tasks)) {
  # task.i = 4
  
  name.task.i <- tasks[task.i]
  
  glms <- unique(glminfo[glminfo$task == name.task.i]$name.glm)
  
  for (glm.i in seq_along(glms)) {
    # glm.i = 1
    
    name.glm.i <- glms[glm.i]
    
    ## get parameter labels

    regressor.labs <- afni(
      "3dinfo", 
      paste0("-label ", 
        file.path(
          dir.analysis, subjs[1], "RESULTS",  name.task.i, paste0(name.glm.i, "_1"),
          paste0("STATS_",  subjs[1], "_1_R_REML.func.gii")
          )
        )
      )
    regressor.labs <- unlist(strsplit(regressor.labs, "\\|"))
    parcs <- names(
      read.delim(
        file.path(
          dir.analysis, subjs[1], "RESULTS",  name.task.i, paste0(name.glm.i, "_1"), 
          paste0("roistats_schaefer400-07_",  subjs[1], "_1_L.txt")
        )
      )[, -(1:2)]
    )

    
    ## build array for task*session

    roimeans <- array(
      NA,
      dim = c(
        length(regressor.labs),
        length(parcs),
        length(subjs),
        2,
        2
      ),
      dimnames = list(term = regressor.labs, parcel = parcs, subj = subjs, run = c("run1", "run2"), hemi = c("L", "R"))
    )

    for (subj.i in seq_along(subjs)) {
      # subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      if (name.subj.i == "432332") next
      
      for (run.i in 1:2) {
        
        for (hemi.i in c("L", "R")) {
          # run.i = 1; hemi.i = "L"
          
          fname.i <- file.path(
            dir.analysis, name.subj.i, "RESULTS",  name.task.i, paste0(name.glm.i, "_", run.i), 
            paste0("roistats_schaefer400-07_",  name.subj.i, "_", run.i, "_", hemi.i, ".txt")
            )
          
          stats.i <- file.path(
            dir.analysis, name.subj.i, "RESULTS",  name.task.i, paste0(name.glm.i, "_", run.i),
            paste0("STATS_",  name.subj.i, "_", run.i, "_", hemi.i, "_REML.func.gii")
          )
          
          labs <- afni("3dinfo", paste0("-label ", stats.i))
          labs <- unlist(strsplit(labs, "\\|"))
          vals <- as.matrix(read.delim(fname.i)[, -(1:2)])
          if (length(labs) != nrow(vals)) stop("bad labels or bad roistats")
          # vals$terms <- labs
          
          is.ok <- 
            identical(labs, regressor.labs) &
            identical(dim(vals), dim(roimeans)[1:2])
          
          if (!is.ok) next
          
          roimeans[, , name.subj.i, run.i, hemi.i] <- vals
          
        }
        
      }
      
    }
   
    saveRDS(roimeans, here("out", "glms", paste0("roistats_", name.task.i, "_", name.glm.i,  ".RDS")))
     
  }
  
}




