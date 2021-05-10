
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
  name.glm = "baseline_aggressive1_EVENTS_censored_shifted",
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)

n.iter <- length(subjs)
pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.iter, clear = FALSE, width = 120
)


M_train <- matrix(1/3, 4, 4)
diag(M_train) <- 0
M_test <- diag(4)


## compute contrast ----

(time.start <- Sys.time())

contrs <- 
  array(
    NA,
    dim = c(n.vert, length(subjs), length(tasks)),
    dimnames = list(vertex = NULL, subj = subjs, task = tasks)
  )
for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 2
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  fname <- here("out", "glms", paste0("sustained_", name.task.i, "_", name.glm.i,  ".RDS"))
  
  ## save if doesn't exist
  if (!file.exists(fname)) {
    
    sustained.i <- read_coefs_2runpm(subjs, name.task.i, name.glm.i, dir.analysis, .getregs = "block#0_Coef")
    saveRDS(sustained.i, fname)
    
  }
  
  betas.i <- readRDS(fname)

  ## remove subj with missing data, remove extra dim:
  betas.i <- betas.i[, 1, !dimnames(betas.i)$subj %in% "432332"]
  contrs[, , name.task.i] <- betas.i
  
}
rm(betas.i)
gc()




## est simil ----


simil <- 
  array(
    NA,
    dim = c(length(tasks), length(parcellation$key), length(subjs)),
    dimnames = list(
      task = tasks, 
      parcel = parcellation$key, 
      subj = subjs
      )
  )


for (subj.i in seq_along(subjs))  {
  # subj.i = 1
  
  res <- enlist(parcellation$key)
  
  
  name.subj.i <- subjs[subj.i]
  contrs.subj.i <- contrs[, subj.i, ]
  
  
  for (parcel.i in seq_along(parcellation$key)) {
    # parcel.i = 20
    
    ## mask:
    
    is.parcel <- schaefer10k == parcel.i
    contrs.subj.parcel.i <- contrs.subj.i[is.parcel, ]
    
    has.signal.all.conds <- !rowSums(abs(contrs.subj.parcel.i) < .Machine$double.eps) > 0
    if (sum(has.signal.all.conds) <= 8) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))
    U <- contrs.subj.parcel.i[has.signal.all.conds, ]  ## get verts with signal
    
    U <- sweep(U, 2, colMeans(U))  ## center
    U_test <- apply(U, "task", function(x) x / sqrt(sum(x^2)))  ## scale patterns to unit length
    U_train <- apply(U_test %*% M_train, 2, function(x) x / sqrt(sum(x^2)))  ## average across folds then re-scale
    
    d <- colSums(U_test * U_train)  ## cosine similarity
    
    simil[, parcel.i, subj.i] <- d
    
  }
  
  pb$tick()  ## progress bar
  
}



## save ----

if (!dir.exists(here("out", "multitask"))) dir.create(here("out", "multitask"))
saveRDS(simil, here("out", "multitask", paste0("taskaxis_aggressive1_cosinesim-unbiased_unpre.RDS")))


(time.run <- Sys.time() - time.start)


