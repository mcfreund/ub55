source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
# suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))



subjs <- subjs[!subjs %in% "432332"]

rsarray <- readRDS(here("out", "multitask", paste0("euclidean-unbiased_unpre.RDS")))


parcels <- dimnames(rsarray)$parcel
n.mods <- length(parcels) * length(subjs)


lab <- mat2vec(rsarray[, , 1, 1])[, c(".row", ".col")]
lab <- lab %>% 
  separate(.row, c("row_task", "row_cond")) %>%
  separate(.col, c("col_task", "col_cond"))
is.bw.task <- lab$row_task != lab$col_task
is.wn.hi <- lab$row_cond == "hi" & lab$col_cond == "hi"
is.wn.lo <- lab$row_cond == "lo" & lab$col_cond == "lo"
is.bw.lo <- lab$row_cond != lab$col_cond
X <- cbind(hihi = is.wn.hi, lolo = is.wn.lo, hilo = is.bw.lo)
X <- X[is.bw.task, ]
X <- sweep(X, 2, colSums(X), "/")



rsa <- function (x) {
  # x = rsarray
  
  ## prepare similarity matrices for fitting ----
  
  ## unwrap into lower-triangle vector
  
  rsvectors <- vector("list", n.mods)
  names(rsvectors) <- combo_paste(subjs, parcels)
  
  for (subj.i in seq_along(subjs)) {
    for (parcel.j in seq_along(parcels)) {
      # subj.i = 1; parcel.j = 1
      
      name.ij <- paste0(subjs[subj.i], "_", parcels[parcel.j])  ## to match name
      v <- mat2vec(x[, , parcel.j, subj.i])$value[is.bw.task]
      rsvectors[[name.ij]] <- v
      
    }
  }
  
  
  # qcor(apply(x, 1:2, mean))
  
  ## fit glms ----
  
  b <- rsvectors %>% map(~ crossprod(., X))
  betas <- as.data.table(do.call(rbind, b))
  betas$id <- names(b)
  betas <- melt(betas, id.vars = "id", value.name = "beta", variable.name = "term")
  
  
  ## format ----
  
  ## create subj, parcel, and hemi cols from id col
  
  stats.subjs <- bind_cols(
    betas,
    reshape2::colsplit(betas$id, pattern = "_", names = c("subj", "parcel"))
  )
  
  
  ## rearrange cols (and drop id col)
  
  stats.subjs %<>% select(subj, parcel, term, beta)
  
}


## run ----

stats.subjs <- rsa(atanh(rsarray))
stats.subjs.univa <- rsa(univa)
stats.subjs.unifo <- rsa(unifo)


## write ----

fwrite(
  stats.subjs,
  here("out", "multitask",  paste0("subjs_jocn.csv"))
)

fwrite(
  stats.subjs.univa,
  here("out", "multitask",  paste0("subjs_jocn_mean.csv"))
)

fwrite(
  stats.subjs.unifo,
  here("out", "multitask",  paste0("subjs_jocn_ssq.csv"))
)

