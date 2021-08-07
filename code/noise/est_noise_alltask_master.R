library(colorout)
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))

glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_aggressive1_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)




for (do.network in c(TRUE, FALSE)) {
  
  if (do.network) {
    rois <- split(parcellation$key, get.network(parcellation$key))
  } else {
    rois <- split(parcellation$key, parcellation$key)
  }
  
  source(here("code", "est_noise_alltask.R"))
  
}