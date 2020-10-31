
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
hi = list(
  Axcpt = "BX",
  Cuedts = c("InConInc", "InConNoInc"),
  Stern = "LL5RN",
  Stroop = c("biasInCon", "PC50InCon")
)
lo = list(
  Axcpt = "BY",
  Cuedts = c("ConInc", "ConNoInc"),
  Stern = "LL5NN",
  Stroop = c("biasCon", "PC50Con")
)


## loop ----

cl <- makeCluster(nrow(glminfo))
registerDoParallel(cl)
time.start <- Sys.time()
z <- foreach(
  glm.i = seq_len(nrow(glminfo)), .inorder = FALSE, .verbose = TRUE,
  .combine = c,
  .packages = c("mikeutils", "here", "data.table", "ggplot2", "grid", "gridExtra", "dplyr")
  ) %dopar% {
  # glm.i = 1
  
  
  res <- enlist(parcellation$key)
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  
  out.dir <- here("out", "fprint", "unimulti_hilo_target_schaefer400-07")
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  
  fig.dir <-file.path(out.dir, "dmat_figs", paste0(name.task.i, "_", name.glm.i))
  if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
  
  ## read betas:
  
  betas.i <- readRDS(
    here::here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )
  
  betas.i <- betas.i[, , , !dimnames(betas.i)$subj %in% "432332", ]  ## remove subj with missing data
  
  ## average across target TRs and calculate contrast:
  
  B_hi <- apply(betas.i[, hi[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean)
  B_lo <- apply(betas.i[, lo[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean)
  B <- B_hi - B_lo
  
  
  
  for (parcel.i in seq_along(parcellation$key)) {
    # parcel.i = 1
    
    is.parcel <- schaefer10k == parcel.i
    B_i <- B[is.parcel, , ]
    
    ## multivariate:
    
    B1 <- t(B_i[, , 1])  ## contrast pattern vectors for each run
    B2 <- t(B_i[, , 2])
    res.mv <- fprint(B1, B2)  ## input dims: subjects*features
    
    
    ##  univariate:
    
    B1_bar <- cbind(rowMeans(B1))  ## get means
    B2_bar <- cbind(rowMeans(B2))
    res.uv <- fprint(B1_bar, B2_bar)
    
    
    ## save 
    
    res[[parcel.i]] <- 
      data.frame(
        multi = res.mv$contrast,
        univa = res.uv$contrast
      )
    
    p <- grid.arrange(
      matplot(res.mv$D) + labs(title = "multivariate"),
      matplot(res.uv$D) + labs(title = "univariate"),
      top = parcellation$key[parcel.i],
      nrow = 1
    )
    
    ggsave(
      file.path(fig.dir, paste0("euclidean_", parcellation$key[parcel.i], ".pdf")), 
      p,
      width = 14, height = 8, units = "cm",
      device = "pdf"
      )
    
     
  }
  
  res <- bind_rows(res, .id = "parcel")
  fwrite(res, file.path(out.dir, paste0("contrast_euclidean_",  name.task.i, "_", name.glm.i, ".csv")))
  
  NULL
  
}
stopCluster(cl)
time.run <- Sys.time() - time.start




