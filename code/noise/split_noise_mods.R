source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
source(here("code", "_atlases.R"))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))
library(mikeutils)

do.network <- TRUE

glminfo <- data.frame(
  task = c("Axcpt", "Cuedts", "Stern", "Stroop"),
  name.glm = c(
    "baseline_aggressive1_EVENTS_censored_shifted"
  ),
  stringsAsFactors = FALSE
)
glminfo <- as.data.table(glminfo)


if (do.network) {
  rois <- unique(get.network(parcellation$key))
} else {
  rois <- parcellation$key
}

fdmasks <- seq(0.3, 0.9, 0.2)



(start.time <- Sys.time())
for (fdmask.i in seq_along(fdmasks)) {
   
   
   ## read inverse covariance matrices:
   
	for (glm.i in seq_len(nrow(glminfo))) {
		# glm.i = 1; fdmask.i = 1; do.runwise = TRUE
    
		name.task.i <- glminfo[glm.i]$task
		name.glm.i <- glminfo[glm.i]$name.glm
    
		d <- readRDS(
		here(
			"out", "invcov", 
			paste0(
				"invcov_", name.task.i, "_", name.glm.i, 
				"_est-runwise", 
				"_parc-", switch(do.network + 1, "parcels400", "network7"), 
				"_cens-fd", fdmasks[fdmask.i], 
				".RDS"
			 )
		  )
		)
	
		for (subj.i in seq_along(subjs)) {
			#subj.i = 1
			
			name.subj.i <- subjs[subj.i]
			
			d_i <- d %>% filter(subj == name.subj.i)
			
			saveRDS(
				d_i, 
				here(
				  "out", "invcov", name.subj.i,
				  paste0(
				   "invcov_", name.task.i, "_", name.glm.i, 
					"_est-runwise", 
					"_parc-", switch(do.network + 1, "parcels400", "network7"), 
					"_cens-fd", fdmasks[fdmask.i], 
					".RDS"
				  )
				)
			)
		
	   }
	   
	}
   
}

(duration <- Sys.time() - start.time)
