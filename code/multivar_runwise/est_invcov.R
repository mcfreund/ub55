
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

options(warn = 2)  ## warning


## input: subjs, task, glmname
## output: RDS file of inverse covariance matrices (array): both hemispheres for subj

dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = "Axcpt",
  name.glm =
    "baseline_Cues_EVENTS_censored_shifted",
  # task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  # name.glm = c(
  #   "baseline_Cues_EVENTS_censored_shifted", 
  #   "baseline_CongruencySwitch_EVENTS_censored_shifted",
  #   "baseline_ListLength_EVENTS_censored_shifted",
  #   "baseline_Congruency_EVENTS_censored_shifted"
  # ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)

errors <- c()


time.start <- Sys.time()
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  X <- readRDS(here("out", "glms", paste0("xmats_", name.task.i, "_", name.glm.i, ".RDS")))
  
  for (subj.i in seq_along(subjs)[1]) {
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
    
    
    
    for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 1
      
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
      
      W_1 <- corpcor::invcov.shrink(E_i1)
      W_2 <- corpcor::invcov.shrink(E_i2)
      
      ## save prewhitening
      
      dir.results.i <- file.path(dir.analysis, name.subj.i, "RESULTS", name.task.i, name.glm.i)
      
      saveRDS(W_1, paste0(dir.results.i, "_1/", "invcov.RDS"))
      saveRDS(W_2, paste0(dir.results.i, "_2/", "invcov.RDS"))
      
      ## save information on data exclusions
      
      saveRDS(
        list(vertex = has.bold, tr = is.included[, 1]), 
        paste0(dir.results.i, "_1/", "inclusions.RDS")
        )
      saveRDS(
        list(vertex = has.bold, tr = is.included[, 2]), 
        paste0(dir.results.i, "_2/", "inclusions.RDS")
      )
      
    }
    
  }
  
}
time.end <- Sys.time() - time.end













  ## 1. loop over subjs, tasks
  ## 2. read resid
  ## 3. loop over ROIs
  ##    - exclude censored frames, vertices with no variance
  ## 4. prewhiten
  ## 5. save (as RDS)



  
  
  # lambda.cov <- attr(W, "lambda")
  # lambda.var <- attr(W, "lambda.var")
  # B %*% (W %*% B)
  
  
  

  
  
  for (hemi.i in c("L", "R")) {
    # hemi.i = "L"
    
    
    
    fname <- file.path(
      .dir, .subjs[subj.i], "RESULTS",  .task, paste0(.glm, "_", run.i),  
      paste0("STATS_", subjs[subj.i], "_", run.i, "_", hemi.i, "_REML.func.gii")
    )
    
    if (!file.exists(fname)) next
    
    B <- mikeutils::read_gifti2matrix(fname)[is.reg, ]
    
    is.ok.i <- isTRUE(all.equal(dim(B), c(n.reg * n.tr, n.vertex)))
    if (!is.ok.i) stop("mismatched beta array")
    
    
    for (reg.i in n.reg) {
      
      is.reg.i <- grepl(regs[reg.i], labs[is.reg])
      B.reg.i <- t(B[is.reg.i, ])
      
      is.ok.ii <- isTRUE(all.equal(dim(betas[inds, reg.i, , subj.i, run.i]), dim(B.reg.i)))
      if (!is.ok.ii) stop("mismatched regressor array")
      
      betas[inds, reg.i, , subj.i, run.i] <- B.reg.i
      
    }
    
  }
      

}









for (glm.i in seq_len(nrow(glminfo))) {
  
  
  
  list.files()
  
  
  
  
  
  
  
  
  
  
  
  # betas.i <- read_betas(subjs, glminfo[glm.i]$task, glminfo[glm.i]$name.glm, dir.analysis )
  # saveRDS(betas.i, here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")))
  
}

