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
  total = nrow(glminfo)*2*5, clear = FALSE, width = 40
)


## functions ----


## TODO:
## add input for FD threshold, read FD, implement
## add input for runwise vs concat (and change code to split by run until end)


invcov <- function(
  name.task.i, name.glm.i, subjs, 
  n.tr, n.vert, n.core,
  rois, parcellation, schaefer10k,
  runwise, fdmask = 0.9
  ) {
  # n.tr = n.trs[name.task.i]
  
  
  ## for each subject (parallel):
  ##   load residual timeseries
  ##   for each roi (serial):
  ##     exclude TRs and vertices
  ##     estimate invcov
  ##     save
  ## output: 
  ##   nested list of invcov matrices (rois within subjects)
  ##     if runwise=TRUE, invcov matrices per run are bound into array
  ##     if runwise=FALSE, one invcov matrix per roi*subject (concatenated along time dimension before estimating)
  
  cl <- makeCluster(n.core)
  registerDoParallel(cl)
  
  res <- foreach(
    subj.i = seq_along(subjs),
    .final = function(x) setNames(x, subjs)
  ) %dopar% {
    # subj.i = 40
    
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
    
    dims.bad <- any(dim(E) != c(n.tr, n.vert))
    if (dims.bad) stop ("bad dims: error time-series")
    
    ## split by run:
    
    E1 <- E[1:(n.tr/2), ]
    E2 <- E[1:(n.tr/2) + n.tr/2, ]
    
    ## load FD timeseries
    
    movregs <- here::here(
      "out", "glms", name.subj.i, "INPUT_DATA", name.task.i, "baseline",
      paste0("Movement_Regressors_", name.task.i, "Bas", c("1_AP.txt", "2_PA.txt"))
    )
    
    if (any(!file.exists(movregs))) return(NA)
    
    movregs1 <- data.table::fread(movregs[1])
    movregs2 <- data.table::fread(movregs[2])
    
    ## get deriv and scale rotations (rads) by 50 (convert to arc length on 50 mm sphere):
    derivs1 <- rbind(0, apply(movregs1, 2, diff)) %*% diag(c(1, 1, 1, 50, 50, 50))
    derivs2 <- rbind(0, apply(movregs2, 2, diff)) %*% diag(c(1, 1, 1, 50, 50, 50))
    
    fd1 <- apply(abs(derivs1), 1, sum)  ## calc FD
    fd2 <- apply(abs(derivs2), 1, sum)
    
    is.lowmot.t1 <- fd1 < fdmask
    is.lowmot.t2 <- fd2 < fdmask
    
    
    ## loop over ROIs
    
    l <- setNames(vector("list", length(rois)), rois)
    for (roi.i in seq_along(rois)) {
      # roi.i = 3
      
      ## mask:
      
      which.parcels <- grep(rois[roi.i], parcellation$key)  ## works with both network and parcel level
      is.roi <- schaefer10k %in% which.parcels
      
      E1_ii <- E1[, is.roi]
      E2_ii <- E2[, is.roi]
      
      ## get vertices that have time-series variance in both runs:
      
      is.good.v1 <- apply(E1_ii, 2, var) > .Machine$double.eps
      is.good.v2 <- apply(E2_ii, 2, var) > .Machine$double.eps
      
      is.good.v <- is.good.v1 & is.good.v2
      
      E1_ii <- E1_ii[, is.good.v]
      E2_ii <- E2_ii[, is.good.v]
      
      ## get trs that have vertex-wise variance, and low-motion in each run (separately):
      
      is.good.t1 <- (apply(E1_ii, 1, var) > .Machine$double.eps) & is.lowmot.t1
      is.good.t2 <- (apply(E2_ii, 1, var) > .Machine$double.eps) & is.lowmot.t2
      
      E1_ii <- E1_ii[is.good.t1, ]
      E2_ii <- E2_ii[is.good.t2, ]
      
      ## estimate invcov
      
      if (runwise) {
        
        
        
        W1_ii <- tryCatch(
          corpcor::invcov.shrink(E1_ii, verbose = FALSE),
          error = NA
          )
        
        W2_ii <- tryCatch(
          corpcor::invcov.shrink(E2_ii, verbose = FALSE),
          error = NA
          )
        
        if (is.na(W1_ii) | is.na(W2_ii)) return(NA)
        
        W_ii <- abind::abind(W1_ii, W2_ii, along = 3)
        
        attr(W_ii, "which.vert") <- which(is.good.v)
        attr(W_ii, "which.tr1") <- which(is.good.t1)
        attr(W_ii, "which.tr2") <- which(is.good.t2)
        for (attr.i in c("lambda", "lambda.estimated", "lambda.var", "lambda.var.estimated")) {
          attr(W_ii, paste0(attr.i, "1")) <- attr(W1_ii, attr.i)
          attr(W_ii, paste0(attr.i, "2")) <- attr(W2_ii, attr.i)
        }
        
        l[[roi.i]] <- W_ii
        
        
      } else {
        
        E_ii <- rbind(E1_ii, E2_ii)
        
        W_ii <- tryCatch(
          corpcor::invcov.shrink(E_ii, verbose = FALSE),
          error = NA
        )
        
        if (is.na(W_ii)) return(NA)
        
        attr(W_ii, "which.vert") <- which(c(is.good.v1, is.good.v2))
        attr(W_ii, "which.tr") <- which(c(is.good.t1, is.good.t2))
        
        l[[roi.i]] <- W_ii
        
      }
      
      
      
      
    }  ## end roi loop
    
    l
    
  }  ## end subj loop
  
  stopCluster(cl)
  
  
  res
  
  
}





## run ----

for (do.runwise in c(TRUE, FALSE)) {
  
  for (fdmask.i in seq(0.3, 0.9, 0.2)) {
  
    for (glm.i in seq_len(nrow(glminfo))) {
      # glm.i = 1; fdmask.i = 0.5; do.runwise = TRUE
      
      
      ## run:
      
      name.task.i <- glminfo[glm.i]$task
      name.glm.i <- glminfo[glm.i]$name.glm
      
      print(name.task.i)
      # (time.start <- Sys.time())
      
      res <- invcov(
        name.task.i = name.task.i, 
        name.glm.i = name.glm.i, 
        subjs = subjs, 
        n.tr = n.trs[name.task.i], 
        n.vert = n.vert, 
        n.core = n.core/2,
        rois = rois, 
        parcellation = parcellation, 
        schaefer10k = schaefer10k,
        runwise = do.runwise, 
        fdmask = fdmask.i
        )
      
      # (time.end <- Sys.time() - time.start)
      pb$tick()  ## progress bar
      
      
      ## wrangle:
      
      d <- tibble(subj = names(res), data = res)  ## to data.frame
      d <- d[lapply(d$data, length) != 1, ]  ## remove bad subjs
      d <- d %>%
        unnest_longer(data, values_to = "invcov", indices_to = "roi") %>%  ## pull out ROIs
        select(subj, roi, invcov)  ## rearrange
      
      ## save:
      
      if (!dir.exists(here("out", "invcov"))) dir.create(here("out", "invcov"))
      
      saveRDS(
        d, 
        here(
          "out", "invcov", 
           paste0(
             "invcov_", name.task.i, "_", name.glm.i, 
             "_est-", switch(do.runwise + 1, "concat", "runwise"), 
             "_parc-", switch(do.network + 1, "parcels400", "network7"), 
             "_cens-fd", fdmask.i, 
             ".RDS"
           )
        )
      )
      
    rm(res, d)
    gc()
      
    }  ## glm
    
  }  ## fdmask
  
}  ## runwise







## test FD thresholds ----



if (FALSE) {  ## NOT RUN
  
  fd <- function(name.task.i, name.glm.i, name.subj.i) {
    
    movregs <- here::here(
      "out", "glms", name.subj.i, "INPUT_DATA", name.task.i, "baseline",
      paste0("Movement_Regressors_", name.task.i, "Bas", c("1_AP.txt", "2_PA.txt"))
    )
    
    if (any(!file.exists(movregs))) return(NA)
    
    movregs1 <- data.table::fread(movregs[1])
    movregs2 <- data.table::fread(movregs[2])
    
    ## get deriv and scale rotations (rads) by 50 (convert to arc length on 50 mm sphere):
    derivs1 <- rbind(0, apply(movregs1, 2, diff)) %*% diag(c(1, 1, 1, 50, 50, 50))
    derivs2 <- rbind(0, apply(movregs2, 2, diff)) %*% diag(c(1, 1, 1, 50, 50, 50))
    
    fd1 <- apply(abs(derivs1), 1, sum)  ## calc FD
    fd2 <- apply(abs(derivs2), 1, sum)
    
    data.table(fd1, fd2)
  
  }
  
  a <- enlist(tasks)
  for (glm.i in seq_len(nrow(glminfo))) {
    # glm.i = 1
    name.task.i <- glminfo[glm.i]$task
    name.glm.i <- glminfo[glm.i]$name.glm
    
    s <- enlist(subjs)
    for (subj.i in seq_along(subjs)) {
      # subj.i = 1
      name.subj.i <- subjs[subj.i]  
      s[[subj.i]] <- fd(name.task.i = name.task.i, name.glm.i = name.glm.i, name.subj.i = name.subj.i)
    }
  
    s[lapply(s, typeof) == "logical"] <- NULL
    a[[glm.i]] <- rbindlist(s, idcol = "subj")
    
  }
  
  
  a <- rbindlist(a, idcol = "task")
  
  fdrates <- a %>%
    pivot_longer(cols = c("fd1", "fd2"), names_to = "run") %>%
    group_by(task, subj, run) %>%
    summarize(
      fd9 = mean(value < 0.9),
      fd8 = mean(value < 0.8),
      fd7 = mean(value < 0.7),
      fd6 = mean(value < 0.6),
      fd5 = mean(value < 0.5),
      fd4 = mean(value < 0.4),
      fd3 = mean(value < 0.3),
      fd2 = mean(value < 0.2),
      fd1 = mean(value < 0.1)
      )
  
  fdrates %>%
    
    melt %>% 
    
    ggplot(aes(variable, value, color = run)) +
    geom_hline(yintercept = 0.5) +
    geom_hline(yintercept = 0.25) +
    geom_point(position = position_dodge(width = 0.2)) +
    facet_grid(rows = vars(task))
  
  
}