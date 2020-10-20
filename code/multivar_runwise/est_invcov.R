
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

# options(warn = 2)  ## warning


## input: subjs, task, glmname
## output: RDS file of inverse covariance matrices (array): both hemispheres for subj

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = "Axcpt",
  name.glm =
    "baseline_Cues_EVENTS_censored_shifted",
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)

errors <- c()

n.iter <- nrow(glminfo) * length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)



time.start <- Sys.time()
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  X <- readRDS(here("out", "glms", paste0("xmats_", name.task.i, "_", name.glm.i, ".RDS")))
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    
    is.included <- X[, "is.included", name.subj.i, ] > 0
    
    E <- abind(
      read_resid(
        .subj = name.subj.i, .task = "Axcpt", .glm = "baseline_Cues_EVENTS_censored_shifted",
        .dir = dir.analysis, .run = 1
        ),
      read_resid(
        .subj = name.subj.i, .task = "Axcpt", .glm = "baseline_Cues_EVENTS_censored_shifted",
        .dir = dir.analysis, .run = 2
        ),
      along = 0
      )
    
    is.tr.ok <- dim(E)[2] == dim(X)[["tr"]]
    if (!is.tr.ok) {
        errors <- c(errors, paste0(name.glm.i, "|", name.subj.i))
        next
    }
    
    
    cl <- makeCluster(n.cores / 2)
    registerDoParallel(cl)
    # time.start <- Sys.time()
    res <- foreach(parcel.i = seq_along(parcellation$key), .inorder = FALSE) %dopar% {
    # for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 20
      
      is.parcel <- schaefer10k == parcel.i
      
      ## mask:
      
      E_i <- E[, , is.parcel]

      ## remove vertices with 0 timecourse variance:
      
      has.bold <- colSums(apply(E_i, c(1, 3), var) > 0) > 1  ## will need to write to file later
      E_i <- E_i[, , has.bold]
      
      ## remove TRs censored from GLM and estimate prewhitening matrix 
      
      E_i1 <- E_i[1, is.included[, 1], ]
      E_i2 <- E_i[2, is.included[, 2], ]
      
      # is.zero <- c(svd(cov(E_i1))$d < 1E-5, svd(cov(E_i2))$d < 1E-5)
      
      W_1 <- tryCatch(
        corpcor::invcov.shrink(E_i1),
        warning = function(x) paste0("warning: ", name.glm.i, "|", name.subj.i, "|", parcel.i),
        error = function(x) paste0("error: ", name.glm.i, "|", name.subj.i, "|", parcel.i)
      )
      W_2 <- tryCatch(
        corpcor::invcov.shrink(E_i2),
        warning = function(x) paste0("warning: ", name.glm.i, "|", name.subj.i, "|", parcel.i),
        error = function(x) paste0("error: ", name.glm.i, "|", name.subj.i, "|", parcel.i)
      )
      
      
      ## save
      
      dir.results.i <- file.path(dir.analysis, name.subj.i, "RESULTS", name.task.i, name.glm.i)
      dir.results.i1 <- paste0(dir.results.i, "_1/", "invcov")
      dir.results.i2 <- paste0(dir.results.i, "_2/", "invcov")
      if (!dir.exists(dir.results.i1)) dir.create(dir.results.i1)
      if (!dir.exists(dir.results.i2)) dir.create(dir.results.i2)
      suffix <- paste0("schaefer400-07_", parcellation$key[parcel.i])
      
      
      had.prob1 <- class(W_1) == "character"
      if (had.prob1) {
        
        data.table::fwrite(
          data.table::as.data.table(W_1), 
          file.path(dir.results.i1, paste0("error_", suffix, ".txt"))
          )
        
      } else {
        
        saveRDS(
          W_1, 
          file.path(dir.results.i1, paste0("invcov_", suffix, ".RDS"))
          )  ## save prewhitening
        saveRDS(
          list(vertex = has.bold, tr = is.included[, 1]), 
          file.path(dir.results.i1, paste0("inclusions_", suffix, ".RDS"))
          )  ## save information on data exclusions
        
      }
      
      
      
      had.prob2 <- class(W_2) == "character" 
      if (had.prob2) {
        
        data.table::fwrite(
          data.table::as.data.table(W_2), 
          file.path(dir.results.i2, paste0("error_", suffix, ".txt"))
          )
        
      } else {
        
        saveRDS(
          W_2, 
          file.path(dir.results.i2, paste0("invcov_", suffix, ".RDS"))
          )
        saveRDS(
          list(vertex = has.bold, tr = is.included[, 2]), 
          file.path(dir.results.i2, paste0("inclusions_", suffix, ".RDS"))
        )
        
      }

    }
    
    
    stopCluster(cl)
    # time.end <- Sys.time() - time.start
    
    pb$tick()  ## progress bar
    
    
  }
  
}
time.run <- Sys.time() - time.start

