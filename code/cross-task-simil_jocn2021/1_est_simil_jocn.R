
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


conditions <- combo_paste(tasks, c("hi", "lo"))

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
  ## remove subj with missing data, get target TRs:
  betas.i <- betas.i[, , target.trs[[glm.i]], !dimnames(betas.i)$subj %in% "432332", ]
  betas.i <- aperm(betas.i, c("vertex", "subj", "reg", "tr", "run"))
  
  betas.i.hi <- rowMeans(betas.i[, , hi[[glm.i]], , ], dims = 2) ## ave across TR, run, relevant regs
  betas.i.lo <- rowMeans(betas.i[, , lo[[glm.i]], , ], dims = 2)
  betas.ii <- abind(betas.i.hi, betas.i.lo, along = 0)
  
  # ## average across target TRs:
  # betas.i <- abind(
  #   apply(betas.i[, hi[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
  #   apply(betas.i[, lo[[glm.i]], target.trs[[glm.i]], , ], c("vertex", "subj", "run"), mean),
  #   along = 0
  # )   ## condition, vertex, subj
  # 
  # betas.ii <- (betas.i[, , , "run1"] + betas.i[, , , "run2"]) / 2  ## average across run
  
  
  ## reshape to match dim(betas):
  
  betas.ii <- aperm(betas.ii, c(2, 1, 3))
  conditions.i <- combo_paste(name.task.i, c("hi", "lo"))  ## run varies slower
  betas[, conditions.i, subjs] <- betas.ii
  
  
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

univa <- 
  array(
    NA,
    dim = c(length(conditions), length(conditions), length(parcellation$key), length(subjs)),
    dimnames = list(.row = conditions, .col = conditions, parcel = parcellation$key, subj = subjs)
  )

unifo <- 
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
    
    univa[, , parcel.i, subj.i] <- tcrossprod(colMeans(B))
    unifo[, , parcel.i, subj.i] <- tcrossprod(colSums(B^2))
    
  }
  
  pb$tick()  ## progress bar
  
  
}



## save ----

if (!dir.exists(here("out", "multitask"))) dir.create(here("out", "multitask"))
saveRDS(simil, here("out", "multitask", paste0("corr-biased_unpre_jocn.RDS")))
saveRDS(univa, here("out", "multitask", paste0("mean-biased_unpre_jocn.RDS")))
saveRDS(unifo, here("out", "multitask", paste0("ssq-biased_unpre_jocn.RDS")))



(time.run <- Sys.time() - time.start)