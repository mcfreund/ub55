source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

glminfo <- data.table(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_null"
  ),
  stringsAsFactors = FALSE
)

resid_type <- "errts"  ## "wherr" or "errts"



## read in design matrices


## get residuals


# subj.i = 1

name.subj.i <- subjs[subj.i]


## read residuals:

E_list <- enlist(tasks)

for (task.i in seq_along(tasks)) {
  # task.i = 1
  
  name.task.i <- glminfo[task.i]$task
  name.glm.i <- glminfo[task.i]$name.glm
  
  for (run.i in c(1, 2)) {
    
    
    
  }
  
  ## load residuals
  
  l <- enlist(c("run1", "run2"))
  for (run.i in 1:2) {  ## test run
    # run.i = 1
    
    eps.name <- here::here(
      "out", "glms", name.subj.i, "RESULTS", name.task.i, paste0(name.glm.i, "_", run.i),
      paste0(resid_type, "_", name.subj.i, "_", run.i, "_", c("L", "R"), "_REML.func.gii")
    )  ## LEFT then RIGHT
    
    if (any(!file.exists(eps.name))) return(NA)
    
    l[[run.i]] <- cbind(
      mikeutils::read_gifti2matrix(eps.name[1]),
      mikeutils::read_gifti2matrix(eps.name[2])
    )  ## LEFT then RIGHT
    
    dims.bad <- any(dim(E_list[[task.i]]) != c(n.trs[name.task.i]/2, n.vert))
    if (dims.bad) stop ("bad dims: error time-series")

  }
  
  
  
  
}
