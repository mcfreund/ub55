## checks events for run-wise GLMs.
## mike freund, 2020-03-31
##



## setup ----


library(here)
library(dplyr)
library(magrittr)
library(data.table)
library(progress)
library(mikeutils)
library(profvis)

## paths

dir.analysis <- "/data/nil-external/ccp/freund/sub-subj-glms/runwise_new"
dir.bluearc.analysis <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/AFNI_ANALYSIS" 

## data

results <- fread(here("freund", "sub-subj-glms", "out", "runwise", paste0("summary_write-events.csv")))

## strings

subjs <- list.files(dir.analysis)
tasks <- c("Axcpt", "Cuedts", "Stern", "Stroop")
sessi <- c("baseline", "proactive", "reactive")

evt.bluearc <- list(
  Axcpt = c("Ang", "AX", "AY", "Bng", "BX", "BY", "button1", "button2", "allTrials"),
  Cuedts = c(
    "ConInc", "ConNoInc", "InConInc", "InConNoInc", "button1", "button2", "allTrials"
  ),
  Stern  = c("LL5NN", "LL5NP", "LL5RN", "not5NN", "not5NP", "not5RN", "button1", "button2", "allTrials"),
  Stroop = c("biasCon", "biasInCon", "PC50Con", "PC50InCon", "buffCon", "allTrials")
)

evt.expected <- list(
  
  Axcpt = 
    c(
      c("Ang", "AX", "AY", "Bng", "BX", "BY") %>% 
      combopaste(c("", "_correct")) %>%
      combopaste(c("_run1", "_run2")),
      combopaste(c("error", "button1", "button2"), c("_run1", "_run2"))
    ),
  
  Cuedts =
    c(
      c("ConInc", "ConNoInc", "InConInc", "InConNoInc") %>% 
      combopaste(c("", "_correct")) %>%
      combopaste(c("_run1", "_run2")),
      combopaste(c("error", "button1", "button2"), c("_run1", "_run2"))
    ),

  Stern = 
    c(
      c("LL5NN", "LL5NP", "LL5RN", "not5NN", "not5NP", "not5RN") %>% 
      combopaste(c("", "_correct")) %>%
      combopaste(c("_run1", "_run2")),
      combopaste(c("error", "button1", "button2"), c("_run1", "_run2"))
    ),
  
  Stroop = 
    c(
      c("biasCon", "biasInCon", "PC50Con", "PC50InCon", "buffCon") %>% 
      combopaste(c("", "_correct")) %>%
      combopaste(c("_run1", "_run2")),
      paste0("error", c("_run1", "_run2"))
    )
  
)

## check block, blockONandOFF, and movreg files ----

allcombs <- expand.grid(subj = subjs, task = tasks, session = sessi)
allcombs$dir <- file.path(dir.analysis, allcombs$subj, "INPUT_DATA", allcombs$task, allcombs$session)

allcombs$blocks.allgood <- FALSE
allcombs$movregs.allgood <- FALSE
allcombs$blocks.exist <- FALSE
allcombs$movregs.exist <- FALSE

for (ii in seq_len(nrow(allcombs))) {
  # ii = 1
  
  comb <- allcombs[ii, ]
  
  fname.block <- 
    combo_paste(
      comb$subj, comb$task, comb$session, 
      c("block", "blockONandOFF"), 
      c("run1.txt", "run2.txt")
    )
  
  fname.movreg <- 
    c(
      paste0("motion_demean_", comb$session, "_", c("run1.1D", "run2.1D")),
      paste0("movregs_FD_mask_", c("run1.txt", "run2.txt"))
    )

  fname.block.abs <- file.path(comb$dir, fname.block)
  fname.movreg.abs <- file.path(comb$dir, fname.movreg)
  
  ## check blocks:
   
  blocks.exist <- all(file.exists(fname.block.abs))
  
  if (blocks.exist) {
    
    allcombs[ii, "blocks.exist"] <- blocks.exist
    
    ## NB: gives warning about incomplete final line... ignore
    blocks <- lapply(fname.block.abs, read.table, header = FALSE, fill = TRUE, na.string = "*", stringsAsFactors = FALSE)
    lens <- vapply(blocks, function(x) length(unlist(x)), FUN.VALUE = numeric(1))
    blocks.allgood <- identical(lens, c(3, 6, 3, 6))
    
    allcombs[ii, "blocks.allgood"] <- blocks.allgood

  }
  
  movregs.exist <- all(file.exists(fname.movreg.abs))
  
  if (movregs.exist) {
    
    allcombs[ii, "movregs.exist"] <- movregs.exist
    
    movregs <- lapply(fname.movreg.abs, read.table, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    ncols <- vapply(movregs, function(x) ncol(x), FUN.VALUE = numeric(1))
    movregs.allgood <- identical(ncols, c(6, 6, 1, 1))
    
    allcombs[ii, "movregs.allgood"] <- movregs.allgood
    
  }
  
}

allcombs %>% filter(!blocks.exist)
allcombs %>% filter(!blocks.allgood)

allcombs %>% filter(!movregs.exist)
allcombs %>% filter(!movregs.allgood)


## build dirs data.frame ----

evt.expected.df <- reshape2::melt(evt.expected, value.name = "evt") %>% rename(task = L1)

dirs <- rbind(
  expand.grid(
    subj = subjs, session = sessi, task = "Axcpt", 
    evt = evt.expected.df %>% filter(task == "Axcpt") %>% pull(evt)
  ),
  expand.grid(
    subj = subjs, session = sessi, task = "Cuedts", 
    evt = evt.expected.df %>% filter(task == "Cuedts") %>% pull(evt)
  ),
  expand.grid(
    subj = subjs, session = sessi, task = "Stern",
    evt = evt.expected.df %>% filter(task == "Stern") %>% pull(evt)
  ),
  expand.grid(
    subj = subjs, session = sessi, task = "Stroop", 
    evt = evt.expected.df %>% filter(task == "Stroop") %>% pull(evt)
  )
)

dirs$rel <- file.path(dirs$subj, "INPUT_DATA", dirs$task, dirs$session)
dirs$abs <- file.path(dir.analysis, dirs$rel)
dirs$fname.rel <- paste0(dirs$subj, "_", dirs$task, "_", dirs$sess, "_", dirs$evt, ".txt")
dirs$fname.abs <- file.path(dirs$abs, dirs$fname.rel)
dirs$run <- ifelse(
  grepl("run1", dirs$evt), "run1",
  ifelse(
    grepl("run2", dirs$evt), "run2",
    ifelse(
      !grepl("run[1-2]", dirs$evt),
      "both", NA
    )
  )
)
dirs$exists <- FALSE
dirs$is.allgood <- FALSE
dirs$is.allna <- FALSE

## remove invalid combinations of regressors (stroop):

invalid <- grepl("buffCon", dirs$evt) & dirs$session %in% c("baseline", "proactive")
dirs <- dirs[!invalid, ]

dirs <- as.data.table(dirs)  ## for speed


## outer loop ----

pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = nrow(dirs), clear = FALSE, width = 120
)


for (subj.i in seq_along(subjs)) {
  # subj.i = which(subjs == "107321")
  
  name.subj.i <- subjs[subj.i]
  is.subj.i <- dirs$subj == name.subj.i
  
  for (task.i in seq_along(tasks)) {
    # task.i = 3
    
    name.task.i <- tasks[task.i]
    is.task.i <- dirs$task == name.task.i
    
    for (sess.i in seq_along(sessi)) {
      # sess.i = 3
      
      name.sess.i <- sessi[sess.i]
      is.sess.i <- dirs$sess == name.sess.i
      
      dirs.i <- dirs[is.sess.i & is.task.i & is.subj.i, ]
      
      
      ## get files from bluearc ----
      
      ## build path
      
      dir.bluearc.i <- file.path(dir.bluearc.analysis, name.subj.i, "INPUT_DATA", name.task.i, name.sess.i)
      fnames.rel.bluearc <- paste0(
        name.subj.i, "_", name.task.i, "_", name.sess.i, "_", evt.bluearc[[name.task.i]], ".txt"
      )
      
      ## add evt info
      ## exception handling with Stroop and buffCon:
      
      if (name.task.i == "Stroop" && name.sess.i %in% c("baseline", "proactive")) {
        
        fnames.rel.bluearc <- fnames.rel.bluearc[-grep("buffCon", fnames.rel.bluearc)]
        evt.bluearc.i <- evt.bluearc[[name.task.i]][-grep("buffCon", evt.bluearc[[name.task.i]])]
        
      } else evt.bluearc.i <- evt.bluearc[[name.task.i]]
      
      fnames.abs.bluearc <- file.path(dir.bluearc.i, fnames.rel.bluearc)
      
      ## read:
      
      # stimtimes.bluearc <- lapply(fnames.abs.bluearc, fread, fill = TRUE)
      stimtimes.bluearc <- lapply(fnames.abs.bluearc, read.table, fill = TRUE, na.strings = "*", header = FALSE)
      names(stimtimes.bluearc) <- evt.bluearc.i
      
      
      ## loop over files ----
      
      for (fname.i in seq_len(nrow(dirs.i))) {
        # fname.i = 2
        # fname.i = which(dirs.i$evt == "not5NN_run2")
        
        name.evt.i <- dirs.i[fname.i]$evt
        name.run.i <- dirs.i[fname.i]$run
        fname.abs.i <- dirs.i[fname.i]$fname.abs
        fname.rel.i <- dirs.i[fname.i]$fname.rel
        
        
        ## read stimtimes ----
        
        is.in.dir <- file.exists(fname.abs.i)
        
        if (is.in.dir) {
          
          dirs.i[fname.i]$exists <- TRUE
          
          # stimtimes <- fread(fname.abs.i, fill = TRUE, header = FALSE, sep = " ")
          stimtimes <- read.table(fname.abs.i, fill = TRUE, na.strings = "*", header = FALSE)
          
          ## is all NA?
          is.allna <- all(unlist(lapply(stimtimes, is.na), use.names = FALSE))
          
          if (is.allna) {
            
            is.allgood <- FALSE
            
          } else {
            
            ## check type
            all.numeric <- all(unlist(lapply(stimtimes, is.numeric), use.names = FALSE))
            if (!all.numeric) stop("not all numeric")
            
            string.evt <- gsub("_correct|_run[1-2]", "", dirs.i[fname.i]$evt)
            ## no "error" files in bluearc, thus read from alltrials:
            if (string.evt == "error") string.evt <- "allTrials"
            
            ## compare to bluearc stimtime:
            
            if (name.run.i %in% c("run1", "run2")) {
              
              index <- switch(name.run.i, run1 = 1, run2 = 2)
              stimtimes.v <- Filter(function(x) !is.na(x), unlist(stimtimes))  ## remove NAs
              
              is.allgood <- all(stimtimes.v %in% stimtimes.bluearc[[string.evt]][index, ])
              
            } else if (name.run.i == "both") {
              
              stimtimes1 <- Filter(function(x) !is.na(x), unlist(stimtimes[1, ]))  ## remove NAs
              stimtimes2 <- Filter(function(x) !is.na(x), unlist(stimtimes[2, ]))  ## remove NAs
              
              is.allgood <- 
                all(stimtimes1 %in% stimtimes.bluearc[[string.evt]][1, ]) & 
                all(stimtimes2 %in% stimtimes.bluearc[[string.evt]][2, ])
              
            }
            
          }
          
          dirs.i[fname.i, c("is.allgood", "is.allna")] <- data.table(is.allgood, is.allna)
          
          rm(is.allgood, is.allna, stimtimes)  ## no accidents
          
        }
        
        rm(is.in.dir)  ## no accidents
        
        pb$tick()  ## progress bar
        # print(fname.i)
        
      }  ## evt loop end
      
      dirs[is.sess.i & is.task.i & is.subj.i, c("exists", "is.allgood", "is.allna")] <-  
        dirs.i[, c("exists", "is.allgood", "is.allna")]
      
      rm(stimtimes.bluearc)  ## no accidents
      
    }  ## sess loop end
    
  }  ## task loop end
  
}  ## subj loop end


fwrite(dirs, here("freund", "sub-subj-glms", "out", "runwise", "summary_stimtime_files_new.csv"))


## examine ----

dirs <- fread(here("freund", "sub-subj-glms", "out", "runwise", "summary_stimtime_files_new.csv"))
results <- fread(here("freund", "sub-subj-glms", "out", "runwise", paste0("summary_write-events.csv")))

sum(!dirs$exists)  ## how many files do not exist?

dirs %>% filter(!exists) %>% select(evt, task) %>% table

## which files have no events?

dirs.allna <- dirs %>% filter(is.allna)
table(dirs.allna$evt, dirs.allna$session)  ## either error evts or _correct evts


## of the files that have events, are there any with !stimtimes %in% stimtimes.bluearc?

sum(!dirs$is.allgood[!dirs$is.allna])
