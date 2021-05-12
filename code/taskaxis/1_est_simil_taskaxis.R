
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


subjs <- subjs[!subjs %in% "432332"]
n.vert <- 20484

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

n.iter <- length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)

## train on separate runs:
# M_train_sub <- matrix(1/3, 4, 4)
# diag(M_train_sub) <- 0
# M_train <- cbind(
#   rbind(M_train_sub, matrix(0, 4, 4)),
#   rbind(matrix(0, 4, 4), M_train_sub),
#   rbind(M_train_sub, matrix(0, 4, 4)),
#   rbind(matrix(0, 4, 4), M_train_sub)
# )
# M_test <- cbind(diag(8), diag(8)[c(5:8, 1:4), ])

## train on both runs:
M_train_sub <- matrix(1/3, 4, 4)
diag(M_train_sub) <- 0
M_train <- cbind(
  rbind(M_train_sub, M_train_sub),  ## tasks vary faster
  rbind(M_train_sub, M_train_sub)
)
M_test <- diag(8)  ## just for consistency with above.






## wrangle betas ----
(time.start <- Sys.time())


contrs <- 
  array(
    NA,
    dim = c(n.vert, length(subjs), length(tasks), 2),
    dimnames = list(vertex = NULL, subj = subjs, task = tasks, run = NULL)
  )
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 2
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  betas.i <- readRDS(
    here::here("out", "glms", paste0("sustained_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )
  
  ## remove subj with missing data, remove extra dim:
  betas.i <- betas.i[, 1, !dimnames(betas.i)$subj %in% "432332", ]
  betas.i <- aperm(betas.i, c("vertex", "subj", "run"))
  
  contrs[, , name.task.i, ] <- betas.i
  
}
rm(betas.i)
gc()




## est simil ----


simil <- 
  array(
    NA,
    dim = c(length(tasks), length(parcellation$key), length(subjs), 2),
    # dim = c(length(tasks), length(parcellation$key), length(subjs), 4),
    dimnames = list(
      task = tasks, 
      parcel = parcellation$key, 
      subj = subjs, 
      # comparison = c("run11", "run22", "run12", "run21")
      comparison = c("run1", "run2")
      )
  )


for (subj.i in seq_along(subjs))  {
  # subj.i = 1
  
  res <- enlist(parcellation$key)
  
  
  name.subj.i <- subjs[subj.i]
  contrs.subj.i <- contrs[, subj.i, , ]
  
  
  for (parcel.i in seq_along(parcellation$key)) {
    # parcel.i = 20
    
    ## mask:
    
    is.parcel <- schaefer10k == parcel.i
    contrs.subj.parcel.i <- contrs.subj.i[is.parcel, , ]
    
    has.signal.all.conds <- !rowSums(abs(contrs.subj.parcel.i) < .Machine$double.eps) > 0
    if (sum(has.signal.all.conds) <= 8) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))
    
    U <- contrs.subj.parcel.i[has.signal.all.conds, , ]  ## get verts with signal
    
    U <- sweep(U, 2:3, colMeans(U))  ## center
    U <- apply(U, c("task", "run"), function(x) x / sqrt(sum(x^2)))  ## scale patterns to unit length
    
    dim(U) <- c(nrow(U), length(tasks)*2)  ## concatenate runs column-wise (tasks vary faster)
    
    U_train <- apply((U %*% M_train), 2, function(x) x / sqrt(sum(x^2)))  ## ave across folds and rescale
    U_test <- (U %*% M_test)
    
    d <- colSums(U_test * U_train)  ## correlate
    
    # dim(d) <- c(4, 4)  ## tasks by comparison (run11, run22, run12, run21)
    dim(d) <- c(4, 2)  ## tasks by comparison (run11, run22, run12, run21)
    
    simil[, parcel.i, subj.i, ] <- d
    
  }
  
  pb$tick()  ## progress bar
  
}



## save ----

if (!dir.exists(here("out", "multitask"))) dir.create(here("out", "multitask"))
# saveRDS(simil, here("out", "multitask", paste0("taskaxis_correlation-unbiased_unpre.RDS")))
saveRDS(simil, here("out", "multitask", paste0("taskaxis_euclidean-unbiased_unpre_train-bothrun.RDS")))



(time.run <- Sys.time() - time.start)
