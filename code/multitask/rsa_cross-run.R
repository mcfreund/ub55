
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))



subjs <- subjs[!subjs %in% "432332"]
conditions <- combo_paste(tasks, c("hi", "lo"), c("run1", "run2"))
n.dim <- length(conditions)
is.lower.tri <- lower.tri(diag(n.dim))
tris <- fread(here("out", "multitask", "rsa-model-tris_cross-run.csv"))
X <- tris %>% select(where(is.numeric)) %>% as.matrix
# X <- scale(X)

rsarray <- readRDS(here("out", "multitask", paste0("corr-biased_unpre.RDS")))



## prepare similarity matrices for regression ----

## check if rows and col names are equal (should be, but just to be sure...)

are.rowcol.equal <- isTRUE(all.equal(dimnames(rsarray)[[1]], dimnames(rsarray)[[2]]))
if(!are.rowcol.equal) stop("rsarray isn't rsarray!")

## get indices and values

parcels <- dimnames(rsarray)$parcel
n.mods <- length(subjs) * length(parcels)

## unwrap into lower-triangle vector

rsvectors <- vector("list", n.mods)
names(rsvectors) <- combo_paste(subjs, parcels)

for (subj.i in seq_along(subjs)) {
  for (parcel.j in seq_along(parcels)) {
    # subj.i = 1; parcel.j = 1
    
    name.ij <- paste0(subjs[subj.i], "_", parcels[parcel.j])  ## to match name
    rsvectors[[name.ij]] <- atanh(rsarray[, , parcel.j, subj.i][is.lower.tri])
    
  }
}

## check numbers

lengths.rsvectors <- map_dbl(rsvectors, length)
if(sum(lengths.rsvectors != 120) > 0) stop("missing row somewhere in rsvectors!")


## fit glms ----

fits <- rsvectors %>% map(~ .lm.fit( x = X, y = .))
betas <- as.data.frame(do.call(rbind, lapply(fits, coef)))

names(betas) <- colnames(X)
betas$id <- rownames(betas)
betas <- melt(as.data.table(betas), id.vars = "id", value.name = "beta", variable.name = "term")


## format ----

## create subj, parcel, and hemi cols from id col

stats.subjs <- bind_cols(
  betas,
  reshape2::colsplit(betas$id, pattern = "_", names = c("subj", "parcel"))
)


## rearrange cols (and drop id col)

stats.subjs %<>% select(subj, parcel, term, beta)

## write ----

fwrite(
  stats.subjs,
  here("out", "multitask",  paste0("subjs_cross-run.csv"))
)

