
R <- array(
  NA,
  c(length(taskruns), length(taskruns), n.subj, length(rois), length(fdmasks), 2, 3),
  list(
    .row = taskruns, .col = taskruns, 
    subj = subjs, roi = rois, 
    fdmask = paste0("mask", fdmasks),
    component = c("var", "cov"),
    measure = c("cor", "cossim", "rmse")
  )
)


(start.time <- Sys.time())
for (fdmask.i in seq_along(fdmasks)) {
  #fdmask.i = 1
  
  for (subj.i in seq_along(subjs)) {
    #subj.i = 1; fdmask.i = 1
    
    name.subj.i <- subjs[subj.i]
    d <- enlist(glminfo$task)
    
    ## read inverse covariance matrices:
    
    for (glm.i in seq_len(nrow(glminfo))) {
      # glm.i = 1; fdmask.i = 1; do.runwise = TRUE
      
      name.task.i <- glminfo[glm.i]$task
      name.glm.i <- glminfo[glm.i]$name.glm
      
      d[[glm.i]] <- readRDS(
        here(
          "out", "invcov", name.subj.i,
          paste0(
            "invcov_", name.task.i, "_", name.glm.i, 
            "_est-runwise", 
            "_parc-", switch(do.network + 1, "parcels400", "network7"), 
            "_cens-fd", fdmasks[fdmask.i], 
            ".RDS"
          )
        )
      )
      
    }
    
    is.bad.task <- vapply(
      d, 
      function(x) {
        x <- c(x$invcov)
        any(is.na(x)) | any(is.null(x)) | length(x) < 1
      },
      logical(1)
    )
    
    if (any(is.bad.task)) next
    
    
    ## loop over ROIs, compute stability stats:
    
    for (roi.i in seq_along(rois)) {
      # roi.i = 1
      
      name.roi.i <- rois[roi.i]
      
      ii_axcpt <- which(d$Axcpt$roi == name.roi.i)
      ii_cuedts <- which(d$Cuedts$roi == name.roi.i)
      ii_stern <- which(d$Stern$roi == name.roi.i)
      ii_stroop <- which(d$Stroop$roi == name.roi.i)
      
      axcpt_vert <- attr(d$Axcpt$invcov[[ii_axcpt]], "which.vert") %>% as.character
      cuedts_vert <- attr(d$Cuedts$invcov[[ii_cuedts]], "which.vert") %>% as.character
      stern_vert <- attr(d$Stern$invcov[[ii_stern]], "which.vert") %>% as.character
      stroop_vert <- attr(d$Stroop$invcov[[ii_stroop]], "which.vert") %>% as.character
      
      good_vert <- Reduce(intersect, list(axcpt_vert, cuedts_vert, stern_vert, stroop_vert)) %>% as.character
      
      if (length(good_vert) < 1) next	
      
      X_cov <- cbind(
        get_lt(d$Axcpt$invcov[[ii_axcpt]], axcpt_vert, good_vert),
        get_lt(d$Cuedts$invcov[[ii_cuedts]], cuedts_vert, good_vert),
        get_lt(d$Stern$invcov[[ii_stern]], stern_vert, good_vert),
        get_lt(d$Stroop$invcov[[ii_stroop]], stroop_vert, good_vert)
      )
      
      X_var <- cbind(
        get_diag(d$Axcpt$invcov[[ii_axcpt]], axcpt_vert, good_vert),
        get_diag(d$Cuedts$invcov[[ii_cuedts]], cuedts_vert, good_vert),
        get_diag(d$Stern$invcov[[ii_stern]], stern_vert, good_vert),
        get_diag(d$Stroop$invcov[[ii_stroop]], stroop_vert, good_vert)
      )
      
      ## get stats:
      
      R[, , name.subj.i, name.roi.i, fdmask.i, "cov", "cor"] <- cor(X_cov)
      R[, , name.subj.i, name.roi.i, fdmask.i, "cov", "cossim"] <- cosinesim_m(X_cov)
      R[, , name.subj.i, name.roi.i, fdmask.i, "cov", "rmse"] <- sqrt(pdist2(t(X_cov), t(X_cov)) / nrow(X_cov))
      
      R[, , name.subj.i, name.roi.i, fdmask.i, "var", "cor"] <- cor(X_var)
      R[, , name.subj.i, name.roi.i, fdmask.i, "var", "cossim"] <- cosinesim_m(X_var)
      R[, , name.subj.i, name.roi.i, fdmask.i, "var", "rmse"] <- sqrt(pdist2(t(X_var), t(X_var)) / nrow(X_var))
      
    }
    
    print(paste0(subj.i, "/", length(subjs), " | ", fdmask.i, "/", length(fdmasks)))
    
  }
  
}
(duration <- Sys.time() - start.time)
saveRDS(R, here("out", "invcov", "invcov-stability_est-runwise.RDS"))  ##Time difference of 6.26064 hours
