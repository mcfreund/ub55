
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

n.iter <- length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)



## wrangle betas ----
(time.start <- Sys.time())


conditions <- combo_paste(tasks, c("hi", "lo"), c("run1", "run2"))

betas <- 
  array(
    NA,
    dim = c(n.vert, length(conditions), length(subjs)),
    dimnames = list(vertex = NULL, condition = conditions, subj = subjs)
  )
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 1
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  ## read betas:
  betas.i <- readRDS(
    here::here("out", "glms", paste0("betas_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS"))
  )
  
  betas.i <- betas.i[, , , !dimnames(betas.i)$subj %in% "432332", ]  ## remove subj with missing data
  
  ## average across target TRs:
  betas.i <- abind(
    apply(betas.i[, hi[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
    apply(betas.i[, lo[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
    along = 0
  )   ## condition, vertex, subj, run
  # names(dimnames(betas.i)) <- c("condition", "vertex", "subj", "run")
  # dimnames(betas.i)$condition <- c("hi", "lo")
  
  ## reshape to match dim(betas):
  
  betas.i <- aperm(betas.i, c(2, 1, 4, 3))  ## vertex, condition, run, subj
  conditions.i <- combo_paste(name.task.i, c("hi", "lo"), c("run1", "run2"))  ## run varies slower
  dim(betas.i) <- c(n.vert, length(conditions.i), length(subjs))  ## vertex, condition_run, subj
  ##.... run varies slower
  # all.equal(betas.i.test[, 1, 1], betas.i[1, , 1, 1])  ## run 1 hi
  # all.equal(betas.i.test[, 2, 1], betas.i[2, , 1, 1])  ## run 1 lo
  # all.equal(betas.i.test[, 3, 1], betas.i[1, , 1, 2])  ## run 2 hi
  # dimnames(betas.i) <- list(vertex = NULL, condition = conditions.i, subj = subjs)
  
  
  betas[, conditions.i, subjs] <- betas.i
  
  
}
rm(betas.i)
gc()





## est simil ----


simil <- 
  array(
    NA,
    dim = c(length(conditions), length(conditions), length(parcellation$key), length(subjs)),
    dimnames = list(.row = conditions, .col = conditions, parcel = parcellation$key, subj = subjs)
  )
for (subj.i in seq_along(subjs))  {
  # subj.i = 1
  
  res <- enlist(parcellation$key)
  
  
  name.subj.i <- subjs[subj.i]
  betas.subj.i <- betas[, , subj.i]
  

  for (parcel.i in seq_along(parcellation$key)) {
    # parcel.i = 20
    
    ## mask:
    
    is.parcel <- schaefer10k == parcel.i
    betas.subj.parcel.i <- betas.subj.i[is.parcel, ]
    
    has.signal.all.conds <- !rowSums(abs(betas.subj.parcel.i) < .Machine$double.eps) > 0
    if (sum(has.signal.all.conds) <= length(conditions)) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))
    
    B <- betas.subj.parcel.i[has.signal.all.conds, ]
    
    
    ## estimate similarity matrices:
    
    simil[, , parcel.i, subj.i] <- cor(B)
    
  }
  
  pb$tick()  ## progress bar

  
}



## save ----

if (!dir.exists(here("out", "multitask"))) dir.create(here("out", "multitask"))
saveRDS(simil, here("out", "multitask", paste0("corr-biased_unpre.RDS")))



(time.run <- Sys.time() - time.start)
