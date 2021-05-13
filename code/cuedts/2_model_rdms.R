source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))



do.network <- TRUE
do.prew <- FALSE
name.glm <- "baseline_cueletnum_EVENTS_censored_shifted"


if (do.network) {
  rois <- unique(get.network(parcellation$key))
} else {
  rois <- parcellation$key
}


simil <- readRDS(
  here("out", "cuedts", 
       paste0(
         "rsamat_cossimil",
         "_est-cval", 
         "_prew-", switch(do.prew + 1, "vanilla", "catwherr"), 
         ".RDS"
       )
  )
)


subjs <- subjs[!subjs %in% c("432332", "DMCC5820265")]
n.subj <- length(subjs)

n.vert <- 20484

regs <- dimnames(simil)[[1]]
n.reg <- length(regs)

n.tr <- dim(simil)[4]

is.lower.tri <- lower.tri(diag(n.reg))

## TODO: add code for reading/writing models

## prepare similarity matrices for regression ----


## get indices and values

rois <- dimnames(simil)$roi
n.mods <- length(subjs) * length(rois) * n.tr

## unwrap into lower-triangle vector

rsvectors <- vector("list", n.mods)
names(rsvectors) <- combo_paste(subjs, rois, paste0("tr", 1:n.tr))

for (subj.i in seq_along(subjs)) {
  for (roi.i in seq_along(rois)) {
    for (tr.i in seq_len(n.tr)) {
      # subj.i = 1; roi.i = 1; tr.i = 1
    
      name.i <- paste0(subjs[subj.i], "_", rois[roi.i], paste0("_tr", tr.i))  ## to match name
      rsvectors[[name.i]] <- scale(simil[, , rois[roi.i], tr.i, subj.i][is.lower.tri])
      
      
    }
  }
}



## fit glms ----


betas <- rsvectors %>% map(~ coef(.lm.fit( x = X, y = .)))
betas <- do.call(rbind, betas)
colnames(betas) <- colnames(X)
betas <- as.data.table(betas, keep.rownames = TRUE)
betas <- rename(betas, id = rn)
betas <- melt(betas, id.vars = "id", value.name = "beta", variable.name = "term")


## format ----

## create subj, parcel, and hemi cols from id col

stats.subjs <- bind_cols(
  betas,
  reshape2::colsplit(betas$id, pattern = "_", names = c("subj", "parcel", "tr"))
)

stats.subjs$tr <- as.numeric(gsub("tr", "", stats.subjs$tr))

## rearrange cols (and drop id col)

stats.subjs %<>% select(subj, parcel, term, tr, beta)

## write ----

fwrite(
  stats.subjs,
  here(
    "out", "cuedts",  
    paste0(
      "stats-subjs",
      "_est-cval", 
      "_prew-", switch(do.prew + 1, "vanilla", "catwherr"), 
      ".RDS"
      )
    )
)



## ---



stats.subjs %>%
  
  # filter(parcel == "Vis") %>%
  
  ggplot(aes(tr, beta)) +
  
  geom_hline(yintercept = 0, color = "black") +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.5) +
  
  facet_grid(rows = vars(parcel), cols = vars(term))


stats.subjs %>%
  
  # filter(parcel == "Vis") %>%
  
  ggplot(aes(tr, beta)) +
  
  geom_hline(yintercept = 0, color = "black") +
  stat_summary(fun = ~t.test(.)$statistic, geom = "line") +
  
  facet_grid(rows = vars(parcel), cols = vars(term))

