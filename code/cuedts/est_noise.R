source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

do.network <- TRUE

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

pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = nrow(glminfo), clear = FALSE, width = 40
)


(time.start <- Sys.time())

for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  name.task.i <- glminfo[glm.i]$task
  name.glm.i <- glminfo[glm.i]$name.glm
  
  print(name.task.i)
  
  ## loop over subjs in parallel ----
  
  cl <- makeCluster(n.core / 2)
  registerDoParallel(cl)
  
  res <- foreach(
    subj.i = seq_along(subjs),
    .final = function(x) setNames(x, subjs)
    ) %dopar% {
    # subj.i = 41
    
    name.subj.i <- subjs[subj.i]
    
    ## load residuals
    
    eps.name <- here::here(
      "out", "glms", name.subj.i, "RESULTS", name.task.i, name.glm.i,
      paste0("wherr_", name.subj.i, "_", c("L", "R"), "_REML.func.gii")
    )  ## LEFT then RIGHT
    
    if (any(!file.exists(eps.name))) return(NA)
    
    E <- cbind(
      mikeutils::read_gifti2matrix(eps.name[1]),
      mikeutils::read_gifti2matrix(eps.name[2])
      )  ## LEFT then RIGHT
    
    dims.bad <- any(dim(E) != c(n.trs[name.task.i], n.vert))
    if (dims.bad) stop ("bad dims: error time-series")
    
    ## loop over ROIs
    
    l <- enlist(rois)
    for (roi.i in seq_along(rois)) {
      # roi.i = 3
      
      ## mask:
      
      which.parcels <- grep(rois[roi.i], parcellation$key)  ## works with both network and parcel level
      is.roi <- schaefer10k %in% which.parcels
      E_ii <- E[, is.roi]
      
      ## get vertices that have time-series variance in both runs:
      
      is.good.v <- apply(E_ii, 2, var) > .Machine$double.eps
      E_ii <- E_ii[, is.good.v]
      
      ## get trs that have vertex-wise variance
      
      is.good.t <- apply(E_ii, 1, var) > .Machine$double.eps
      E_ii <- E_ii[is.good.t, ]
      
      ## estimate invcov
      
      W_ii <- corpcor::invcov.shrink(E_ii, verbose = FALSE)
      attr(W_ii, "which.vert") <- which(is.good.v)
      attr(W_ii, "which.tr") <- which(is.good.t)
      
      
      l[[roi.i]] <- W_ii
      
    }  ## end roi loop
    
    l
  
  }  ## end subj loop
  
  stopCluster(cl)
  
  pb$tick()  ## progress bar
  
  
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
           "invcov_", name.task.i, "_", name.glm.i, 
           "_est-concat", 
           "_parc-", switch(do.network + 1, "parcels400", "network7"), 
           ".RDS"
         )
      )
    )
    
  }
  
  
  
}


(time.stop <- Sys.time() - time.start)
