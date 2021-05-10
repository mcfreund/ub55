
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


do.network <- TRUE  ## metwork or parcel level?


# subjs <- subjs[!subjs %in% "432332"]
n.vert <- 20484

glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm.block = c(
    "baseline_Cues_EVENTS_censored_shifted",
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
  ),
  name.glm.noblock = c(
    "baseline_Cues_EVENTS_censored_shifted_noblock",
    "baseline_CongruencySwitch_EVENTS_censored_shifted_noblock",
    "baseline_ListLength_EVENTS_censored_shifted_noblock",
    "baseline_Congruency_EVENTS_censored_shifted_noblock"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)

if (do.network) {
  rois <- unique(get.network(parcellation$key))
} else {
  rois <- parcellation$key
}


n.iter <- length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)


## wrangle betas ----


contrs <- 
  array(
    NA,
    dim = c(n.vert, length(subjs), length(tasks), 2),
    dimnames = list(vertex = NULL, subj = subjs, task = tasks, run = NULL)
  )
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 2
  
  name.glm.i <- glminfo[glm.i]$name.glm.block
  name.task.i <- glminfo[glm.i]$task
  
  betas.i <- readRDS(
    here::here("out", "glms", paste0("sustained_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm.block,  ".RDS"))
  )
  
  betas.i <- betas.i[, 1, , ]  ## remove extra dim
  betas.i <- aperm(betas.i, c("vertex", "subj", "run"))
  
  contrs[, , name.task.i, ] <- betas.i
  
}
rm(betas.i)
gc()




## project ----


if (!dir.exists(here("out", "taskaxis"))) dir.create(here("out", "taskaxis"))

(time.start <- Sys.time())


for (task.i in seq_along(tasks)) {
  # task.i = 4
  
  name.task.i <- glminfo[task.i]$task
  contrs.task.i <- contrs[, , name.task.i, ]
  
  ## get n.tr:
  
  n.tr <- n.trs[name.task.i] / 2
  
  proj <-
    array(
      NA,
      dim = c(n.tr, 4, length(rois), length(subjs)),
      dimnames = list(
        tr = NULL,
        fold = c("1_1", "2_1", "1_2", "2_2"),  ## train, test
        roi = rois,
        subj = subjs
      )
    )
  
  
  for (subj.i in seq_along(subjs))  {
    # subj.i = 1
    
    res <- enlist(parcellation$key)
    
    if (subjs[subj.i] == "DMCC5820265") next
    
    name.subj.i <- subjs[subj.i]
    U <- contrs.task.i[, name.subj.i, ]

    
    for (run.i in 1:2) {  ## test run
      # run.i = 1
      
      eps.name <- here::here(
        "out", "glms", name.subj.i, "RESULTS", name.task.i, paste0(glminfo[task.i]$name.glm.noblock, "_", run.i),
        paste0("errts_", name.subj.i, "_", run.i, "_", c("L", "R"), "_REML.func.gii")
      )  ## LEFT then RIGHT
      eps <- cbind(read_gifti2matrix(eps.name[1]), read_gifti2matrix(eps.name[2]))  ## LEFT then RIGHT
      
      for (mask.i in seq_along(rois)) {
        # mask.i = 5
        
        ## mask:
        
        which.parcels <- grep(rois[mask.i], parcellation$key)
        is.roi <- schaefer10k %in% which.parcels
        U_roi <- U[is.roi, ]  ## task axis
        x <- eps[, is.roi]  ## TR-level patterns
        
        ## exclude verts with no signal:
        
        has.signal.all.runs <- !rowSums(abs(U_roi) < .Machine$double.eps) > 0
        good.vert <- apply(x, 2, var) > .Machine$double.eps
        j <- good.vert & has.signal.all.runs  ## subset by
        
        x_j <- x[, j]
        U_roi_j <- U_roi[j, ]
        
        # good.tr <- apply(x, 1, var) > .Machine$double.eps  ## TRs that weren't censored

        ## project:
        
        U_roi_j <- apply(U_roi_j, 2, function(x) x / sqrt(sum(x^2)))  ## scale vectors to unit length
        p_j <- x_j %*% U_roi_j
        
        
        proj[, paste0(1:2, "_", run.i), mask.i, subj.i] <- p_j
        
        
      }
      
    }
    
    
  }
  
  
  saveRDS(
    proj, 
    here(
      "out", "taskaxis", 
      paste0("projections_task-", name.task.i, "_parc-network_prew-vanilla_glm-noblock_resi-errts.RDS")
      )
    )
  
  print(paste0(name.task.i, " done."))
  
}

(time.run <- Sys.time() - time.start)












## ----


# 
# mu <- apply(proj, 1:3, mean, na.rm = TRUE)
# mu <- apply(mu[, c(2, 3), ], c(1, 3), mean, na.rm = TRUE)
# 
# 
# data.frame(p = mu[, "Vis"], tr = 1:540) %>%
#   
#   ggplot(aes(x = tr, y = p)) +
#   geom_line()
# 

# 
# dir.subj <- here("out", "glms", subjs)
# d <- data.table(expand.grid(subj = subjs, task = tasks))
# d$n.missing <- as.numeric(NA)
# for (subj.i in seq_along(dir.subj)) {
#   
#   for (task.i in seq_along(tasks)) {
#     
#     dirs <- list.files(file.path(dir.subj[subj.i], "RESULTS", tasks[task.i]), pattern = "noblock", full.names = TRUE)
#     
#     d[subj == subjs[subj.i] & task == tasks[task.i]]$n.missing <- 4 - sum(grepl("wherr", list.files(dirs)))
#     
#   }
#   
# }
# 
# d[n.missing > 0]
# 
