## about ----
##
## saves xmats as arrays within RDS files.



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
source(here("code", "_vars.R"))


## build data.frame for looping

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = c("Stroop"),
  name.glm = c(
    "fix-LSA_EVENTS_censored"
    )
)
dirs <- full_join(dirs, glminfo, by = "task")
dirs$name.glm <- paste0(dirs$session, "_", dirs$name.glm)
dirs$run1 <- file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, paste0(dirs$name.glm, "_1"), "X_1.xmat.1D")
dirs$run2 <- file.path(dir.analysis, dirs$subj, "RESULTS", dirs$task, paste0(dirs$name.glm, "_2"), "X_2.xmat.1D")
dirs <- pivot_longer(dirs, cols = c("run1", "run2"), names_to = "run", values_to = "fname.xmat")
dirs$exists.xmat <- file.exists(dirs$fname.xmat)

dirs$fname.censor <- dirs$fname.xmat %>% gsub("RESULTS", "INPUT_DATA", .) %>% gsub("/baseline_.*", "/baseline/movregs_FD_mask_", .)
dirs$fname.censor <- paste0(dirs$fname.censor, dirs$run, ".txt")

dirs <- filter(dirs, !subj %in% "432332")  ## for now, remove those that don't exist

dirs <- as.data.table(dirs)  ## for fast extracting



## loop ----


## initialize list

# xmats <- setNames(vector("list", 4), tasks)

pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = nrow(glminfo)*length(unique(dirs$subj)), clear = FALSE, width = 120
)

for (task.i in seq_along(unique(glminfo$task))) {
  # task.i = 4
  
  name.task.i <- as.character(unique(glminfo$task)[task.i])
  
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
        tr = n.trs[[name.task.i]] / 2, 
        regressor = length(xlabels) + 1, 
        subj = length(subjs),
        run = 2
      ),
      dimnames = list(tr = NULL, regressor = c(xlabels, "is.included"), subj = subjs, run = c("run1", "run2"))
    )
    
    for (subj.i in seq_along(subjs)) {
      # subj.i = 1
      
      name.subj.i <- subjs[subj.i]
      
      if (name.subj.i == "432332") next
      
      fname.i <- dirs[subj == name.subj.i & task == name.task.i & name.glm == name.glm.i]$fname.xmat
      
      xmat.subj.i1 <- read_xmat(fname.i[1])
      xmat.subj.i2 <- read_xmat(fname.i[2])
      
      fname.censor.i <- dirs[subj == name.subj.i & task == name.task.i & name.glm == name.glm.i]$fname.censor
      
      xmat.subj.i1 <- cbind(xmat.subj.i1, as.matrix(fread(fname.censor.i[1])))
      xmat.subj.i2 <- cbind(xmat.subj.i2, as.matrix(fread(fname.censor.i[2])))
      colnames(xmat.subj.i1)[ncol(xmat.subj.i1)] <- "is.included"
      colnames(xmat.subj.i2)[ncol(xmat.subj.i2)] <- "is.included"

      ## check for match
      
      are.matching.xmats <- 
        identical(dimnames(xmat.subj.i1), dimnames(xmat.subj.i2)) && 
        identical(dim(xmat.subj.i1), dim(xmat.subj.i2))
      
      if (are.matching.xmats) {
    
        are.matched.xlabels <- identical(colnames(xmat.i), colnames(xmat.subj.i1))
        are.matched.dims <- all(dim(xmat.i)[c("tr", "regressor")] == dim(xmat.subj.i1))
        
        if (!(are.matched.xlabels && are.matched.dims)) stop("unexpected xlabels or dims")
        
      } else stop("mismatched xlabels or dims across runs")
      
      
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

