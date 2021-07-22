
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
} else {
  rois <- parcellation$key
}


B <- readRDS(here("out", "glms", paste0("betas_Cuedts_", name.glm,  ".RDS")))


subjs <- subjs[!subjs %in% c("432332", "DMCC5820265")]
n.subj <- length(subjs)

n.vert <- 20484

regs <- dimnames(B)$reg
n.reg <- length(regs)

n.tr <- dim(B)[3]

M <- contrast_matrix(n.reg, regs)


pb <- progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.subj, clear = FALSE, width = 120
)





## est simil ----



simil <- 
  array(
    NA,
    dim = c(n.reg, n.reg, length(rois), n.tr, length(subjs)),
    dimnames = list(regs, regs, roi = rois, tr = NULL, subj = subjs)
  )



(time.start <- Sys.time())


for (subj.i in seq_along(subjs))  {
  # subj.i = 1
  
  
  name.subj.i <- subjs[subj.i]
  B_i <- B[, , , name.subj.i, ]
  
  
  ## load residuals
  
  if (do.prew) {
    
    E <- array(NA, dim = c(n.trs["Cuedts"]/2, n.vert, 2))
    
    for (run.i in 1:2) {
      # run.i = 1
      
      eps.name <- here(
        "out", "glms", name.subj.i, "RESULTS", "Cuedts", paste0(name.glm, "_", run.i),
        paste0("wherr_", name.subj.i, "_", run.i, "_", c("L", "R"), "_REML.func.gii")
      )  ## LEFT then RIGHT
      
      # inds <- (run.i-1)*n.trs["Cuedts"]/2 + 1:(n.trs["Cuedts"]/2)  ## lec
      
      E[, , run.i] <- cbind(read_gifti2matrix(eps.name[1]), read_gifti2matrix(eps.name[2]))  ## LEFT then RIGHT
      
    }
    
  }
  
  ## loop over ROIs

  for (roi.i in seq_along(rois)) {
    # roi.i = 3
    
    ## mask:
    
    which.parcels <- grep(rois[roi.i], parcellation$key)  ## works with both network and parcel level
    is.roi <- schaefer10k %in% which.parcels
    
    B_ii <- B_i[is.roi, , , ]
    
    
    ## remove bad vertices
    
    is.good.v <- !rowSums(abs(B_ii) < .Machine$double.eps) > 0
    # if (sum(has.signal.all.conds) <= n.roi) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))

    ## estimate invcov
    
    
    if (do.prew) {
      
      E_ii <- E[, is.roi, ]
      
      ## get vertices that have time-series variance in both runs:
      
      has.e.all.runs <- rowSums(apply(E_ii, 2:3, var) > .Machine$double.eps) == 2
      is.good.v <- is.good.v & has.e.all.runs
      
      E_ii <- E_ii[, is.good.v, ]
      
      ## get trs that have vertex-wise variance
      
      E_ii <- rbind(E_ii[, , 1], E_ii[, , 2])  ## concatenate runs
      
      is.good.t <- apply(E_ii, 1, var) > .Machine$double.eps
      
      E_ii <- E_ii[is.good.t, ]
      
      ## estimate invcov
      
      W_ii <- corpcor::invcov.shrink(E_ii)
      
    }
    
    
    U <- B_ii[is.good.v, , , ]
    
    
    ## scale
    
    U_unit <- sweep(
      U, 2:4, apply(U, c("reg", "tr", "run"), function(x) sqrt(sum(x^2))), "/"
    )  ## to unit length
    
    
    for (tr.i in seq_len(n.tr)) {
      # tr.i = 1
      
      U_i <- U_unit[, , tr.i, ]
      
      
      if (do.prew) {
        
        simil[, , rois[roi.i], tr.i, name.subj.i] <- cvdist(U_i, M, W = W_ii)
        
      } else {
        
        simil[, , rois[roi.i], tr.i, name.subj.i] <- cvdist(U_i, M)
        
      }
      
      
    }  ## end tr loop
    
    
  }  ## end roi loop
  
  pb$tick()  ## progress bar
  
}  ## end subj loop




## save ----

if (!dir.exists(here("out", "cuedts"))) dir.create(here("out", "cuedts"))
saveRDS(
  simil, 
  here("out", "cuedts", 
       paste0(
         "rsamat_cossimil",
         "_est-cval", 
         "_parc-", switch(do.network + 1, "parcels400", "network7"), 
         "_prew-", switch(do.prew + 1, "vanilla", "catwherr"), 
         ".RDS"
         )
       )
  )


(time.run <- Sys.time() - time.start)


