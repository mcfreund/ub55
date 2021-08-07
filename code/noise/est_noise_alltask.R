

## loop over subjs in parallel ----

cl <- makeCluster(n.core / 2)
registerDoParallel(cl)

res <- foreach(
  subj.i = seq_along(subjs),
  .final = function(x) setNames(x, subjs),
  .verbose = TRUE,
  .packages = c("data.table", "here", "mikeutils", "purrr", "corpcor")
) %dopar% {
  # subj.i = 1
  
  name.subj.i <- subjs[subj.i]
  
  
  ## read residuals:
  
  E_list <- enlist(tasks)
  for (task.i in seq_along(tasks)) {
    # task.i = 1
    
    name.task.i <- glminfo[task.i]$task
    name.glm.i <- glminfo[task.i]$name.glm
    
    ## load residuals
    
    eps.name <- here(
      "out", "glms", name.subj.i, "RESULTS", name.task.i, name.glm.i,
      paste0("wherr_", name.subj.i, "_", c("L", "R"), "_REML.func.gii")
    )  ## LEFT then RIGHT
    
    
    if (any(!file.exists(eps.name))) return(NA)
    
    E_list[[task.i]] <- cbind(
      read_gifti2matrix(eps.name[1]),
      read_gifti2matrix(eps.name[2])
    )  ## LEFT then RIGHT

    dims.bad <- any(dim(E_list[[task.i]]) != c(n.trs[name.task.i], n.vert))
    if (dims.bad) stop ("bad dims: error time-series")
    
  }
  
  
  ## loop over ROIs
  
  l <- enlist(names(rois))
  for (roi.i in seq_along(l)) {
    # roi.i = 3
    
    ## mask:
    name.roi.i <- names(rois)[roi.i]
    
    which.parcels <- match(rois[[name.roi.i]], parcellation$key)  ## ACTUALLY works with both network and parcel level
    ### which.parcels <- grep(rois[[name.roi.i]], parcellation$key)  ## BUG! do not un-comment! 2021-08-07
    
    is.roi <- schaefer10k %in% which.parcels
    E_ii <- map(E_list, ~.x[, is.roi])
    
    ## get vertices that have time-series variance in both runs:
    ## NB: filter by those that have time-series variance in EACH run? --- may need to split runs
    is.good.v <- rowSums(map_df(E_ii, ~apply(.x, 2, var) > .Machine$double.eps)) == 4
    
    E_ii <- map(E_ii, ~.x[, is.good.v])
    
    ## get trs that have vertex-wise variance
    is.good.t <- map(E_ii, ~apply(.x, 1, var) > .Machine$double.eps)
    E_ii <- map2(E_ii, is.good.t, ~ .x[.y, ])
    
    ## scale and concatenate
    
    # E_ii <- map(E_ii, scale)  ## don't scale, variances relatively consistent across tasks (see stabil. analysis)
    E_iic <- do.call(rbind, E_ii)

    ## estimate invcov
    
    W_ii <- invcov.shrink(E_iic, verbose = FALSE)
    attr(W_ii, "which.vert") <- which(is.good.v)
    attr(W_ii, "which.tr") <- is.good.t
    
    
    l[[roi.i]] <- W_ii
    
  }  ## end roi loop
  
  l

} ## end subj loop
  
stopCluster(cl)
(time.stop <- Sys.time() - time.start)



## wrangle:

d <- 
  tibble(subj = names(res), data = res) %>%
  unnest_longer(data, values_to = "invcov") %>%  ## pull out ROIs
  select(subj, roi = data_id, invcov)  ## rearrange


## save:
  
for (name.subj.i in subjs) {
  
  if (!dir.exists(here("out", "invcov", name.subj.i))) dir.create(here("out", "invcov", name.subj.i))
  
  saveRDS(
    d %>% filter(subj == name.subj.i), 
    here(
      "out", "invcov", name.subj.i,
      paste0(
        "invcov_alltask_baseline_aggressive1_EVENTS_censored_shifted", 
        "_est-concat", 
        "_parc-", switch(do.network + 1, "parcels400", "network7"), 
        ".RDS"
      )
    )
  )
  
}


