
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))
library(mikeutils)

do.network <- TRUE
do.calcstabil <- FALSE

glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_aggressive1_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)


if (do.network) {
  rois <- unique(get.network(parcellation$key))
} else {
  rois <- parcellation$key
}


## functions ----

pdist2 <- function(A,B) {
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  tmp - 2 * tcrossprod(A,B)
}


cosinesim_m <- function(M) {
  S <- diag(1 / apply(M, 2, function(x) sqrt(sum(x^2))))
  MS <- tcrossprod(M, S)
  crossprod(MS)
}


get_lt <- function(m, vert_names, vert_names_keep, do.runwise = TRUE) {
  
  if (do.runwise) {
    
    dimnames(m) <- list(row = vert_names, col = vert_names, run = c("run1", "run2"))
    
    m <- m[vert_names_keep, vert_names_keep, ]
    
    m1 <- m[, , 1]
    m2 <- m[, , 2]
    
    m1 <- m1[lower.tri(m1)]
    m2 <- m2[lower.tri(m2)]
    
    cbind(m1, m2)
    
  } else {
    
    dimnames(m) <- list(row = vert_names, col = vert_names)
    
    m <- m[vert_names_keep, vert_names_keep]
    
    m <- m[, ]
    
    m[lower.tri(m)]
    
  }
  
  
}

get_diag <- function(m, vert_names, vert_names_keep, do.runwise = TRUE) {
  
  if (do.runwise) {
    
    dimnames(m) <- list(row = vert_names, col = vert_names, run = c("run1", "run2"))
    
    m <- m[vert_names_keep, vert_names_keep, ]
    
    m1 <- m[, , 1]
    m2 <- m[, , 2]
    
    m1 <- diag(m1)
    m2 <- diag(m2)
    
    cbind(m1, m2)
    
  } else {
    
    dimnames(m) <- list(row = vert_names, col = vert_names)
    
    m <- m[vert_names_keep, vert_names_keep]
    
    diag(m)
    
  }
  
  
}




## runwise noise estimates ----
R <- array(
  NA,
  c(length(tasks), length(tasks), n.subj, length(rois), 2, 3),
  list(
    .row = tasks, .col = tasks, 
    subj = subjs, roi = rois, 
    component = c("var", "cov"),
    measure = c("cor", "cossim", "rmse")
  )
)


(start.time <- Sys.time())
for (subj.i in seq_along(subjs)) {
  #subj.i = 1
  
  name.subj.i <- subjs[subj.i]
  d <- enlist(glminfo$task)
  
  ## read inverse covariance matrices:
  
  for (glm.i in seq_len(nrow(glminfo))) {
    # glm.i = 1
    
    name.task.i <- glminfo[glm.i]$task
    name.glm.i <- glminfo[glm.i]$name.glm
    
    d[[glm.i]] <- readRDS(
      here(
        "out", "invcov", name.subj.i,
        paste0(
          "invcov_", name.task.i, "_", name.glm.i, 
          "_est-concat", 
          "_parc-", switch(do.network + 1, "parcels400", "network7"), 
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
      get_lt(d$Axcpt$invcov[[ii_axcpt]], axcpt_vert, good_vert, do.runwise = FALSE),
      get_lt(d$Cuedts$invcov[[ii_cuedts]], cuedts_vert, good_vert, do.runwise = FALSE),
      get_lt(d$Stern$invcov[[ii_stern]], stern_vert, good_vert, do.runwise = FALSE),
      get_lt(d$Stroop$invcov[[ii_stroop]], stroop_vert, good_vert, do.runwise = FALSE)
    )
    
    X_var <- cbind(
      get_diag(d$Axcpt$invcov[[ii_axcpt]], axcpt_vert, good_vert, do.runwise = FALSE),
      get_diag(d$Cuedts$invcov[[ii_cuedts]], cuedts_vert, good_vert, do.runwise = FALSE),
      get_diag(d$Stern$invcov[[ii_stern]], stern_vert, good_vert, do.runwise = FALSE),
      get_diag(d$Stroop$invcov[[ii_stroop]], stroop_vert, good_vert, do.runwise = FALSE)
    )
    
    ## get stats:
    
    R[, , name.subj.i, name.roi.i, "cov", "cor"] <- cor(X_cov)
    R[, , name.subj.i, name.roi.i, "cov", "cossim"] <- cosinesim_m(X_cov)
    R[, , name.subj.i, name.roi.i, "cov", "rmse"] <- sqrt(pdist2(t(X_cov), t(X_cov)) / nrow(X_cov))
    
    R[, , name.subj.i, name.roi.i, "var", "cor"] <- cor(X_var)
    R[, , name.subj.i, name.roi.i, "var", "cossim"] <- cosinesim_m(X_var)
    R[, , name.subj.i, name.roi.i, "var", "rmse"] <- sqrt(pdist2(t(X_var), t(X_var)) / nrow(X_var))
    
  }
  
  print(paste0(subj.i, "/", length(subjs)))
  
}
(duration <- Sys.time() - start.time)
saveRDS(R, here("out", "invcov", "invcov-stability_est-concat.RDS"))  ##Time difference of 6.26064 hours
