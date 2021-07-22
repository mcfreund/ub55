
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



subj.i <- 54
glm.i = 1
tr.i <- 12
name.reg.i <- "biasInCon"
name.reg.j <- "biasCon"
parcel.i <- 379

simil <- readRDS(
  here(
    "out", "rsa", paste0("euclidean-cv-unsta-unpre_", glminfo[glm.i]$task, "_", glminfo[glm.i]$name.glm,  ".RDS")
  )
)
name.glm.i <- glminfo[glm.i]$name.glm
name.task.i <- glminfo[glm.i]$task
fname.L1 <- file.path(
  dir.analysis, subjs[subj.i], "RESULTS",  name.task.i, paste0(name.glm.i, "_", 1),  
  paste0("STATS_", subjs[subj.i], "_", 1, "_L_REML.func.gii")
)
fname.R1 <- file.path(
  dir.analysis, subjs[subj.i], "RESULTS",  name.task.i, paste0(name.glm.i, "_", 1),  
  paste0("STATS_", subjs[subj.i], "_", 1, "_R_REML.func.gii")
)

fname.L2 <- file.path(
  dir.analysis, subjs[subj.i], "RESULTS",  name.task.i, paste0(name.glm.i, "_", 2),  
  paste0("STATS_", subjs[subj.i], "_", 2, "_L_REML.func.gii")
)
fname.R2 <- file.path(
  dir.analysis, subjs[subj.i], "RESULTS",  name.task.i, paste0(name.glm.i, "_", 2),  
  paste0("STATS_", subjs[subj.i], "_", 2, "_R_REML.func.gii")
)


BL1 <- read_gifti2matrix(fname.L1)
BR1 <- read_gifti2matrix(fname.R1)
B1 <- cbind(BL1, BR1)

BL2 <- read_gifti2matrix(fname.L2)
BR2 <- read_gifti2matrix(fname.R2)
B2 <- cbind(BL2, BR2)


labels <- afni("3dinfo", paste0("-label ", fname.L1))
labels <- strsplit(labels, "|", fixed = TRUE)[[1]]
this.reg.i <- which(labels == paste0(name.reg.i, "#", tr.i - 1, "_Coef"))
this.reg.j <- which(labels == paste0(name.reg.j, "#", tr.i - 1, "_Coef"))
is.parcel <- schaefer10k == parcel.i



dir.results <- file.path(dir.analysis, subjs[subj.i], "RESULTS", name.task.i, name.glm.i)
suffix <- paste0("schaefer400-07_", parcellation$key[parcel.i])
W <- list(run1 = NA, run2 = NA)
inclusions <- list(run1 = NA, run2 = NA)
for (run.i in 1:2) {
  # run.i = 1
  
  dir.results.run <- paste0(dir.results, "_", run.i, "/", "invcov")
  W[[run.i]] <- readRDS(file.path(dir.results.run, paste0("invcov_", suffix, ".RDS")))
  inclusions[[run.i]] <- readRDS(file.path(dir.results.run, paste0("inclusions_", suffix, ".RDS")))
  
}



manual <- (B1[this.reg.i, is.parcel][inclusions$run1$vertex] - B1[this.reg.j, is.parcel][inclusions$run1$vertex]) %*%
(B2[this.reg.i, is.parcel][inclusions$run1$vertex] - B2[this.reg.j, is.parcel][inclusions$run1$vertex]) / sum(inclusions$run1$vertex)
manual <- c(manual)


# simil <- simil[, , , , subjs[!subjs %in% "448347"]]
auto <- simil[name.reg.i, name.reg.j, tr.i, parcel.i, subj.i]
dimnames(simil[, , tr.i, parcel.i, ])$subj[subj.i]

all.equal(manual, auto)

all.equal(
  c(B2[this.reg.i, is.parcel][inclusions$run1$vertex] %*% pracma::sqrtm(W$run2)$B), 
  c(pracma::sqrtm(W$run2)$B %*%  B2[this.reg.i, is.parcel][inclusions$run1$vertex])
  )

plot(
  B2[this.reg.i, is.parcel][inclusions$run1$vertex],
  pracma::sqrtm(W$run2)$B %*% B2[this.reg.i, is.parcel][inclusions$run1$vertex]
)


plot(
  pracma::sqrtm(W$run1)$B[lower.tri(pracma::sqrtm(W$run1)$B, diag = TRUE)],
  pracma::sqrtm(W$run2)$B[lower.tri(pracma::sqrtm(W$run2)$B, diag = TRUE)]
)

plot(
  pracma::sqrtm(W$run1)$B[lower.tri(pracma::sqrtm(W$run1)$B)],
  pracma::sqrtm(W$run2)$B[lower.tri(pracma::sqrtm(W$run2)$B)]
)


W1 <- pracma::sqrtm(W$run1)$B
W2 <- pracma::sqrtm(W$run2)$B


W$run2 %>%
  
  symmat4ggplot %>%
  
  ggplot(aes(v1, v2, fill = value)) +
  geom_raster() +
  
  scale_fill_viridis_c()



solve(W$run1) %>%
   
  cov2cor %>%
  
  symmat4ggplot %>%
  
  ggplot(aes(v1, v2, fill = value)) +
  geom_raster() +
  
  scale_fill_viridis_c()

plot(svd(solve(W$run2))$d)

plot(diag(solve(W$run2)), diag(solve(W$run1)))


