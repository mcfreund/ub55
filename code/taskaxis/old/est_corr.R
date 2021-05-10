source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

subjs <- subjs[!subjs %in% "432332"]


dirs <- expand.grid(subj = subjs, task = tasks, session = "baseline", stringsAsFactors = FALSE)
glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),#, "Stroop"),
  name.glm = c(
    "baseline_Cues_EVENTS_censored_shifted", 
    "baseline_CongruencySwitch_EVENTS_censored_shifted",
    "baseline_ListLength_EVENTS_censored_shifted",
    "baseline_Congruency_EVENTS_censored_shifted"
    # "baseline_fix-item_EVENTS_censored"
  )
)
glminfo <- as.data.table(glminfo)

networks <- unique(get.network(parcellation$key))

errors <- c()


for (glm.i in seq_len(nrow(glminfo))) {
  # glm.i = 4
  
  name.glm.i <- glminfo[glm.i]$name.glm
  name.task.i <- glminfo[glm.i]$task
  
  
  X <- readRDS(here("out", "glms", paste0("xmats_", name.task.i, "_", name.glm.i, ".RDS")))
  
  
  S <- array(
    NA,
    dim = c(
      tr_run1 = dim(X)[["tr"]],
      tr_run2 = dim(X)[["tr"]],
      parcel  = length(networks), 
      subj    = length(subjs)
    ),
    dimnames = list(
      tr_run1 = NULL,
      tr_run1 = NULL,
      network  = networks,
      subj    = subjs
    )
  )
  
  
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1
    
    name.subj.i <- subjs[subj.i]
    
    is.included <- X[, "is.included", name.subj.i, ] > 0
    
    E <- abind(
      read_resid(
        .subj = name.subj.i, .task = name.task.i, .glm = name.glm.i,
        .dir = dir.analysis, .run = 1
      ),
      read_resid(
        .subj = name.subj.i, .task = name.task.i, .glm = name.glm.i,
        .dir = dir.analysis, .run = 2
      ),
      along = 0
    )
    
    is.tr.ok <- dim(E)[2] == dim(X)[["tr"]]
    if (!is.tr.ok) {
      errors <- c(errors, paste0(name.glm.i, "|", name.subj.i))
      next
    }
    
    
    for (parcel.i in seq_along(networks)) {
      # parcel.i = 4
      
      is.parcel <- schaefer10k %in% which(networks[parcel.i] == get.network(parcellation$key))
      
      ## mask:
      
      E_i <- E[, , is.parcel]
      
      ## remove vertices with 0 timecourse variance:
      
      has.bold <- colSums(apply(E_i, c(1, 3), var) > 0) > 1  ## will need to write to file later
      E_i <- E_i[, , has.bold]

      S[, , parcel.i, subj.i] <- cor(t(E_i[1, , ]), t(E_i[2, , ]))
      
      
    }
    
  }
  
  saveRDS(
    S,
    here(
      "out", "taskcoding", paste0("correlation-crossrun_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
    )
  )
  
  print(glm.i)
  
}








