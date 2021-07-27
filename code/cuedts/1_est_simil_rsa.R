#!/usr/bin/env Rscript


## argparser ----
# https://bitbucket.org/djhshih/argparser/src/master/

library(argparser, quietly = TRUE)


# Create a parser
p <- arg_parser("Estimate cross-validated distances")


# Add command line arguments

p <- add_argument(p, "--prew", help = "Should patterns be spatially prewhitened?", flag = TRUE)
p <- add_argument(p, "--network", help = "Should network (versus parcel) level be used?", flag = TRUE)


p <- add_argument(
  p, "--name_out", 
  help = "Name of .RDS output file ('' if prew=resid; 'Task-shortglmname-prewname' if prew=invcov, e.g., Axcpt-aggressive1-catwherr", 
  type = "character"
  )

p <- add_argument(
  p, "--name_glm",
  help = "Name of first-level GLM (not including task, e.g., baseline_cueletnum_EVENTS_censored_shifted)", 
  type = "character", 
  default = "baseline_cueletnum_EVENTS_censored_shifted"
  )

p <- add_argument(
  p, "--name_invcov",
  help = "Name of invcov .RDS file to read (in out/invcov/./.). 
  If incov should be estimated from same first-level GLM as betas, input 'estimate'.", 
  type = "character"
)



## Parse the command line arguments
argv <- parse_args(p)


if (interactive()) {
  
  argv$prew <- TRUE
  argv$name_out <- "Axcpt_aggressive1-catwherr"
  argv$name_invcov <- "invcov_Axcpt_baseline_aggressive1_EVENTS_censored_shifted_est-concat"
  argv$network <- TRUE
  
}



## input validation:

if (argv$prew == FALSE & !is.na(argv$name_invcov)) stop("--name_incov specified but --prew not flagged.")



## setup ----

source(here::here("code", "_vars.R"))
library(cifti)
library(magrittr)
library(tibble)
source(here::here("code", "_atlases.R"))
source(here::here("code", "_funs.R"))


## get rois:

if (argv$network) {
  rois <- unique(get.network(parcellation$key))
} else {
  rois <- parcellation$key
}


<<<<<<< HEAD
subjs <- subjs[!subjs %in% c("432332", "DMCC5820265")]  ## exclude some subjs
=======
B <- readRDS(here("out", "glms", paste0("betas_Cuedts_", name.glm,  ".RDS")))
>>>>>>> a254b73a1259023929990fb6a4cd247026a30ca0


## load data and get info from it

B <- readRDS(here::here("out", "glms", paste0("betas_Cuedts_", argv$name_glm,  ".RDS")))

regs <- dimnames(B)$reg
n.reg <- length(regs)
n.tr <- dim(B)[3]


M <- mikeutils::contrast_matrix(n.reg, regs)  ## contrast matrix for RSA


pb <- progress::progress_bar$new(
  format = " running [:bar] :percent eta: :eta (elapsed: :elapsed)",
  total = n.subj, clear = FALSE, width = 120
)


cat(unlist(argv), "\n")



<<<<<<< HEAD

=======
>>>>>>> a254b73a1259023929990fb6a4cd247026a30ca0
## est simil ----



simil <-
  array(
    NA,
    dim = c(n.reg, n.reg, length(rois), n.tr, length(subjs)),
    dimnames = list(regs, regs, roi = rois, tr = NULL, subj = subjs)
  )



time.start <- Sys.time()
# cat(time.start, "\n")




for (subj.i in seq_along(subjs))  {
  # subj.i = 1


  name.subj.i <- subjs[subj.i]
  B_i <- B[, , , name.subj.i, ]


  if (argv$prew) {

    ## (load residuals & estimate noise model) | (load noise model):

    if (!is.na(argv$name_invcov) & argv$name_invcov == "estimate") {

      E <- array(NA, dim = c(n.trs["Cuedts"]/2, n.vert, 2))

      for (run.i in 1:2) {
        # run.i = 1

        eps.name <- here::here(
          "out", "glms", name.subj.i, "RESULTS", "Cuedts", paste0(argv$name_glm, "_", run.i),
          paste0("wherr_", name.subj.i, "_", run.i, "_", c("L", "R"), "_REML.func.gii")
        )  ## LEFT then RIGHT

        # inds <- (run.i-1)*n.trs["Cuedts"]/2 + 1:(n.trs["Cuedts"]/2)  ## lec

        E[, , run.i] <- cbind(
          mikeutils::read_gifti2matrix(eps.name[1]),
          mikeutils::read_gifti2matrix(eps.name[2])
          )  ## LEFT then RIGHT

      }

    } else {

      W_l <- readRDS(
        here::here(
          "out", "invcov", subjs[subj.i],
          paste0(argv$name_invcov, "_parc-", switch(argv$network + 1, "parcels400", "network7"), ".RDS")
          )
      )

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

    is.good.v.beta <- !rowSums(abs(B_ii) < .Machine$double.eps) > 0
    # if (sum(has.signal.all.conds) <= n.roi) stop(paste0("bad parcel: ", subj.i, " ", parcel.i))

    ## read / estimate invcov

    if (argv$prew) {

      if (argv$name_invcov == "estimate") {

        E_ii <- E[, is.roi, ]

        ## get vertices that have time-series variance in both runs:

        has.e.all.runs <- rowSums(apply(E_ii, 2:3, var) > .Machine$double.eps) == 2
        is.good.v <- is.good.v.beta & has.e.all.runs
        E_ii <- E_ii[, is.good.v, ]

        ## get trs that have vertex-wise variance

        E_ii <- rbind(E_ii[, , 1], E_ii[, , 2])  ## concatenate runs
        is.good.t <- apply(E_ii, 1, var) > .Machine$double.eps
        E_ii <- E_ii[is.good.t, ]

        ## estimate invcov

        W_ii <- corpcor::invcov.shrink(E_ii)

      } else {

        W <- W_l$invcov[[which(W_l$roi == rois[roi.i])]]
        is.good.v.W <- is.good.v.beta[attr(W, "which.vert")]  ## of length nrow(W)
        W_ii <- W[is.good.v.W, is.good.v.W]

        is.good.v.invcov <- 1:length(is.good.v.beta) %in% attr(W, "which.vert")
        is.good.v <- is.good.v.invcov & is.good.v.beta  ## of length n_vert in ROI

      }

    }


    U <- B_ii[is.good.v, , , ]  ## if prew, will additionally subset by prew is.good.verts


    ## scale

    U_unit <- sweep(
      U, 2:4, apply(U, c("reg", "tr", "run"), function(x) sqrt(sum(x^2))), "/"
    )  ## to unit length


    for (tr.i in seq_len(n.tr)) {
      # tr.i = 1

      U_i <- U_unit[, , tr.i, ]

      if (argv$prew) {
        simil[, , rois[roi.i], tr.i, name.subj.i] <- mikeutils::cvdist(U_i, M, W = W_ii)
      } else {
        simil[, , rois[roi.i], tr.i, name.subj.i] <- mikeutils::cvdist(U_i, M)
      }

    }



  }  ## end roi loop



  pb$tick()  ## progress bar

}  ## end subj loop





## save ----

if (!dir.exists(here::here("out", "cuedts"))) dir.create(here::here("out", "cuedts"))
saveRDS(
  simil,
  here::here("out", "cuedts",
       paste0(
         "rsamat_cossimil",
         "_est-cval",
         "_parc-", switch(argv$network + 1, "parcels400", "network7"),
         "_prew-", switch(argv$prew + 1, "vanilla", argv$name_out),
         ".RDS"
         )
       )
  )



time.run <- Sys.time() - time.start


cat(time.run, "\n")

