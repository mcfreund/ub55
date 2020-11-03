
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

subjs <- subjs[!subjs %in% "432332"]


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



subj.i <- 23
glm.i = 4
run.i <- 1 
tr.i <- 1
name.reg.i <- "PC50Con"
hemi.i <- "L"



name.glm.i <- glminfo[glm.i]$name.glm
name.task.i <- glminfo[glm.i]$task


## read betas:

betas.i <- readRDS(
  here::here("out", "glms", paste0("betas_", name.task.i, "_", name.glm.i,  ".RDS"))
)

fname <- file.path(
  dir.analysis, subjs[subj.i], "RESULTS",  name.task.i, paste0(name.glm.i, "_", run.i),  
  paste0("STATS_", subjs[subj.i], "_", run.i, "_", hemi.i, "_REML.func.gii")
)
B <- read_gifti2matrix(fname)
labels <- afni("3dinfo", paste0("-label ", fname))
labels <- strsplit(labels, "|", fixed = TRUE)[[1]]
name.reg <- paste0(name.reg.i, "#", tr.i - 1, "_Coef")
this.reg <- which(labels == name.reg)

# dim(B)
head(B[this.reg, ])
# dim(betas.i)
betas.i[1:6, name.reg.i, tr.i, subjs[subj.i], run.i]

all.equal(
  betas.i[1:10242, name.reg.i, tr.i, subjs[subj.i], run.i],
  B[this.reg, ]
)
 
# plot(
#   betas.i[1:10242, name.reg.i, tr.i, subjs[subj.i], run.i],
#   B[this.reg, ]
# )
# 
# 


## checking ROIstats.
# 
# d <- setNames(vector("list", nrow(glminfo)), glminfo$task)
# for (glm.i in seq_len(nrow(glminfo))) {
#   # glm.i = 1
#   
#   a <- readRDS(here("out", "glms", paste0("roistats_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm, ".RDS")))
#   a <- a[, , subjs[!subjs %in% "432332"], ]
#   a <- as.data.table(reshape2::melt(a))
#   
#   a <- separate(a, term, c("term", "stat"), "_")
#   a <- a[!grepl("Full|block", term) & stat != "Fstat"]
#   a <- separate(a, term, c("term", "tr"), "#")
#   
#   a$network <- a$parcel %>% gsub("LH_|RH_", "", .) %>% gsub("_.*", "", .)
#   a$hemi <- substr(a$parcel, 1, 1)
#   
#   a <- dcast(a, term + tr + parcel + subj + run + hemi + network ~ stat, value.var = "value")
#   names(a) <- tolower(names(a))
#   
#   d[[glm.i]] <- a
#   
# }
# rm(glm.i, a)
# 
# ## fix cuedts labels: incorrect! (mislabeled in afni code)
# d$Cuedts$term <- gsub("NoInc$", "Repeat", d$Cuedts$term)
# d$Cuedts$term <- gsub("Inc$", "Switch", d$Cuedts$term)
# 
# dim(betas.i)
# dim(a)
# 
# betas.i[, ]
# 
# 



