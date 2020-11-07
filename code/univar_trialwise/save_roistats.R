
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))


## input: subjs, task, session, glmname
## output: RDS files

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = c("Stroop"),
  name.glm = c(
    "baseline_fix-LSA_EVENTS_censored"
  )
)
glminfo <- as.data.table(glminfo)

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
    
    
    ## build array for task*session
    
    roimeans <- array(
      NA,
      dim = c(
        length(regressor.labs),
        length(parcellation$key),
        length(subjs),
        2
      ),
      dimnames = list(term = regressor.labs, parcel = parcellation$key, subj = subjs, run = c("run1", "run2"))
    )
    
    for (subj.i in seq_along(subjs)) {
      # subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      if (name.subj.i == "432332") next
      
      for (run.i in 1:2) {
        
        for (hemi.i in c("L", "R")) {
          # run.i = 1; hemi.i = "L"
          
          if (hemi.i == "L") {
            inds <- 1:200
          } else inds <- 201:400
          
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
          if (!identical(labs, regressor.labs)) stop("bad regressor labels")
          if (200 != ncol(vals)) stop("bad regressor labels")
          
          roimeans[, inds, name.subj.i, run.i] <- vals
          
        }
        
      }
      
    }
    
    saveRDS(roimeans, here("out", "glms", paste0("roistats_", name.task.i, "_", name.glm.i,  ".RDS")))
    
  }
  
}
