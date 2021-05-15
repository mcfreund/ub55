source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))



do.network <- FALSE
do.prew <- TRUE
name.glm <- "baseline_cueletnum_EVENTS_censored_shifted"


if (do.network) {
  rois <- unique(get.network(parcellation$key))
  fname.parc <- ""
} else {
  rois <- parcellation$key
  fname.parc <- "_parc-parcels400"
}


simil <- readRDS(
  here("out", "cuedts", 
       paste0(
         "rsamat_cossimil" ,
         "_est-cval", 
         fname.parc,
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


ctsmods <- readRDS(here::here("out", "cuedts", "rsmods_full.RDS"))
X <- as.matrix(data.table::fread(here::here("out", "cuedts", "rsmods_lt.csv")))


## fit models ----


## prep:

vecs <- apply(simil, 3:5, function(x) scale(x[is.lower.tri]))  ## vectorize and scale
vecs <- as.data.table(vecs)


## regress:

lmfit <- function(x, y, nms = colnames(X)) {
  
  b <- coef(.lm.fit(x = X, y = y))
  
  data.frame(b = b, term = nms)
  
}

stats.subjs <- vecs %>%
  group_by(roi, tr, subj) %>%
  nest %>%
  mutate(data = map(data, ~lmfit(X, .$value))) %>%
  unnest(cols = data)


## write ----

fwrite(
  stats.subjs,
  here(
    "out", "cuedts",  
    paste0(
      "stats-subjs",
      "_est-cval",
      fname.parc,
      "_prew-", switch(do.prew + 1, "vanilla", "catwherr"), 
      ".RDS"
      )
    )
)

