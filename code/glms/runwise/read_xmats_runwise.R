## about ----



## setup ----


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
source(here("code", "_strings.R"))


## build data.frame for looping

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop", "Stroop"),
  name.glm = c(
    "Cues_EVENTS_censored_shifted", 
    "CongruencyIncentive_EVENTS_censored_shifted", 
    "ListLength_EVENTS_censored_shifted",
    "Congruency_EVENTS_censored_shifted",
    "fix-item_EVENTS_censored_shifted"
    )
)
dirs <- full_join(dirs, glminfo, by = "task")
dirs$name.glm <- paste0(dirs$session, "_", dirs$name.glm)
dirs$run1 <- file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, paste0(dirs$name.glm, "_1"), "X_1.xmat.1D")
dirs$run2 <- file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, paste0(dirs$name.glm, "_2"), "X_2.xmat.1D")
dirs <- pivot_longer(dirs, cols = c("run1", "run2"), names_to = "run", values_to = "fname.xmat")
dirs$exists.xmat <- file.exists(dirs$fname.xmat)

dirs <- filter(dirs, exists.xmat)  ## for now, remove those that don't exist

dirs <- as.data.table(dirs)  ## for fast extracting


## loop ----


## initialize list

# xmats <- setNames(vector("list", 4), tasks)

pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = nrow(dirs), clear = FALSE, width = 120
)

for (task.i in seq_along(tasks)) {
  # task.i = 4
  
  name.task.i <- tasks[task.i]
  
  glms <- unique(dirs[task == name.task.i]$name.glm)
  
  for (glm.i in seq_along(glms)) {
    # glm.i = 1
    
    name.glm.i <- glms[glm.i]
    
    ## get parameter labels
    
    sample.fname <- dirs[task == name.task.i & name.glm == name.glm.i]$fname.xmat[1]
    xlabels <- afni("1d_tool.py", paste0("-infile ", sample.fname, " -show_group_labels"))
    xlabels <- gsub("(.*) label (.*)", "\\2", xlabels)

    ## build array for task*session
    
    xmat.i <- array(
      NA,
      dim = c(
        tr = n.trs[[name.task.i]], 
        regressor = length(xlabels), 
        subj = length(subjs),
        run = 2
      ),
      dimnames = list(tr = NULL, regressor = xlabels, subj = subjs, run = c("run1", "run2"))
    )
    
    for (subj.i in seq_along(subjs)) {
      # subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      fname.i <- dirs[subj == name.subj.i & task == name.task.i]$fname.xmat
      
      xmat.subj.i1 <- read_xmat(fname.i[1])
      xmat.subj.i2 <- read_xmat(fname.i[2])
      
      ## check for match
      
      if (identical(dimnames(xmat.subj.i1), dimnames(xmat.subj.i2)) && identical(dim(xmat.subj.i1), dim(xmat.subj.i2))) {
    
        are.matched.xlabels <- identical(colnames(xmat.i), colnames(xmat.subj.i1))
        are.matched.dims <- all(dim(xmat.i)[c("tr", "regressor")] == dim(xmat.subj.i1) * c(2, 1))
        if (!(are.matched.xlabels && are.matched.dims)) {
        #   xmat.subj.i1 <- as.numeric(NA)
        #   xmat.subj.i2 <- as.numeric(NA)
        #   print(noquote(paste0("skipping subj ", name.subj.i)))
          stop("mismatched")
        }
        
      } else stop("mismatched")
      
      
      xmat.i[, , name.subj.i, 1] <- xmat.subj.i1
      xmat.i[, , name.subj.i, 2] <- xmat.subj.i2
      
      rm(xmat.subj.i1, xmat.subj.i2)  ## no accidents
      
      pb$tick()  ## progress bar
      
    }
    
    saveRDS(xmat.i, here("out", "glms", paste0("xmats_",name.task.i, "_", name.glm.i,  ".RDS")))
    
    rm(xmat.i)
    gc()
    
  }

}

