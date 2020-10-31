
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

subjs <- subjs[!subjs %in% "432332"]


glminfo <- data.frame(
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



## loop ----

# n.iter <- nrow(glminfo) * length(subjs)
# pb <- progress_bar$new(
#   format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
#   total = n.iter, clear = FALSE, width = 120
# )



cl <- makeCluster(nrow(glminfo))
registerDoParallel(cl)
time.start <- Sys.time()
res <- foreach(
  glm.i = seq_len(nrow(glminfo)), .inorder = FALSE, .verbose = TRUE,
  .combine = c,
  .packages = c("mikeutils", "here", "data.table")
  ) %dopar% {
# for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  
  ## read betas:
  
  betas.i <- readRDS(
    here::here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )

  ## make array for similarity matrices:

  regressors <- dimnames(betas.i)$reg
  trs <- 1:dim(betas.i)[which(names(dimnames(betas.i)) == "tr")]
  
  simil <- array(
    NA,
    dim = c(
      .row   = length(regressors),
      .col   = length(regressors),
      tr     = length(trs),
      parcel = length(parcellation$key), 
      subj   = length(subjs),
      prewh  = 2,
      stand  = 2
    ),
    dimnames = list(
      .row   = regressors,
      .col   = regressors,
      tr     = trs,
      parcel = parcellation$key, 
      subj   = subjs,
      prewh  = c("true", "false"),
      stand  = c("true", "false")
    )
  )
  
  cmat <- mikeutils::contrast_matrix(length(regressors), regressors)  ## contrast matrix
  
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    betas.subj.i <- betas.i[, , , subj.i, ]
    
    ## loop over parcels:
    
    for (parcel.i in seq_along(parcellation$key)) {
      # parcel.i = 20
      
      dir.results <- file.path(dir.analysis, name.subj.i, "RESULTS", name.task.i, name.glm.i)
      suffix <- paste0("schaefer400-07_", parcellation$key[parcel.i])
      
      
      ## read whitening matrices:
      
      W <- list(run1 = NA, run2 = NA)
      inclusions <- list(run1 = NA, run2 = NA)
      for (run.i in 1:2) {
        # run.i = 1
        
        dir.results.run <- paste0(dir.results, "_", run.i, "/", "invcov")
        W[[run.i]] <- readRDS(file.path(dir.results.run, paste0("invcov_", suffix, ".RDS")))
        inclusions[[run.i]] <- readRDS(file.path(dir.results.run, paste0("inclusions_", suffix, ".RDS")))

      }
      
      W1 <- pracma::sqrtm(W$run1)$B  ## root
      W2 <- pracma::sqrtm(W$run2)$B
      
      
      ## mask:
      
      is.parcel <- schaefer10k == parcel.i
      betas.subj.parcel.i <- betas.subj.i[is.parcel, , , ]
      
      inclusions.are.ok <- 
        identical(inclusions$run1$vertex, inclusions$run2$vertex) & 
        isTRUE(all.equal(length(inclusions$run1$vertex), dim(betas.subj.parcel.i)[1]))
      if (!inclusions.are.ok) stop("inclusions not ok")
      
      betas.subj.parcel.i <- betas.subj.parcel.i[inclusions$run1$vertex, , , ]  ## exclude verts with 0 var(BOLD)
      
      B <- aperm(betas.subj.parcel.i, c(2, 1, 3, 4))  ## condition*vertex*run
      
      ## estimate similarity matrices:
      
      for (tr.i in trs) {
        # tr.i = 1
        
        B1 <- B[, , tr.i, "run1"]
        B2 <- B[, , tr.i, "run2"]
        
        D             <- distance_cv(B1, B2, cmat, regressors)
        D_prewh       <- distance_cv(B1 %*% W1, B2 %*% W2, cmat, regressors)
        D_stand       <- distance_cv(t(scale(t(B1))), t(scale(t(B2))), cmat, regressors)
        D_stand_prewh <- distance_cv(t(scale(t(B1))) %*% W1, t(scale(t(B2))) %*% W2, cmat, regressors)
        
        simil[, , tr.i, parcel.i, subj.i, 2, 2] <- D
        simil[, , tr.i, parcel.i, subj.i, 1, 2] <- D_prewh
        simil[, , tr.i, parcel.i, subj.i, 2, 1] <- D_stand
        simil[, , tr.i, parcel.i, subj.i, 1, 1] <- D_stand_prewh
        
      }

    }
    
    rm(D, D_prewh, D_stand, D_stand_prewh, W, W1, W2, B1, B2, B)

    # pb$tick()  ## progress bar
    
    
  }
  
  
  ## save
  
  if (!dir.exists(here("out", "rsa"))) dir.create(here("out", "rsa"))
  # saveRDS(simil, here("out", "rsa", paste0("simil_cv_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")))
  
  saveRDS(
    simil[, , , , , "true", "true"],
    here(
      "out", "rsa", paste0("euclidean-cv-stand-prewh_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
    )
  )
  
  saveRDS(
    simil[, , , , , "false", "false"],
    here(
      "out", "rsa", paste0("euclidean-cv-unsta-unpre_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
    )
  )
  
  
  saveRDS(
    simil[, , , , , "true", "false"],
    here(
      "out", "rsa", paste0("euclidean-cv-stand-unpre_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
    )
  )
  
  saveRDS(
    simil[, , , , , "false", "true"],
    here(
      "out", "rsa", paste0("euclidean-cv-unsta-prewh_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
    )
  )
  
  
  NULL
  
}
stopCluster(cl)
time.end <- Sys.time() - time.start




# for (glm.i in seq_len(nrow(glminfo))) {
#   
#   simil <- readRDS(
#     here("out", "rsa", paste0("simil_cv_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
#     )
# 
#   saveRDS(
#     simil[, , , , , "true", "true"],
#     here(
#       "out", "rsa", paste0("euclidean-cv-stand-prewh_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
#     )
#   )
#   
#   saveRDS(
#     simil[, , , , , "false", "false"],
#     here(
#       "out", "rsa", paste0("euclidean-cv-unsta-unpre_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
#     )
#   )
# 
#   
#   saveRDS(
#     simil[, , , , , "true", "false"],
#     here(
#       "out", "rsa", paste0("euclidean-cv-stand-unpre_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
#     )
#   )
#   
#   saveRDS(
#     simil[, , , , , "false", "true"],
#     here(
#       "out", "rsa", paste0("euclidean-cv-unsta-prewh_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
#     )
#   )
# 
# }
