## make models -----

regs <- c(
  "let_odd_con", "let_odd_vow", "let_eve_con", "let_eve_vow", "num_odd_con",
  "num_odd_vow", "num_eve_con", "num_eve_vow"
)

x_cue <- model.matrix(~as.factor(substr(regs, 1, 3)) + 0)
x_num <- model.matrix(~as.factor(substr(regs, 5, 7)) + 0)
x_let <- model.matrix(~as.factor(substr(regs, 9, 11)) + 0)
x_stim <- model.matrix(~as.factor(substr(regs, 5, 11)) + 0)

x_resp <- cbind(
  l = grepl("let.*vow|num_eve", regs),
  r = grepl("let.*con|num_odd", regs)
)

x_tt <- cbind(
  incon = grepl("eve_con|odd_vow", regs),
  congr = grepl("eve_vow|odd_con", regs)
)


cue <- tcrossprod(x_cue, x_cue)
num <- tcrossprod(x_num, x_num)
let <- tcrossprod(x_let, x_let)
stim <- tcrossprod(x_stim, x_stim)

cuelet <- tcrossprod(x_cue[, 1], x_cue[, 1])
cuenum <- tcrossprod(x_cue[, 2], x_cue[, 2])
tar <- cuelet * let + cuenum * num  ## target stim coding
dis <- cuelet * num + cuenum * let  ## distractor stim coding

resp <- tcrossprod(x_resp, x_resp)

incon <- tcrossprod(x_tt[, "incon"], x_tt[, "incon"])
congr <- tcrossprod(x_tt[, "congr"], x_tt[, "congr"])
iicc <- tcrossprod(x_tt, x_tt)

X <- cbind(
  cue = cue[lower.tri(cue)],
  # stim = stim[lower.tri(cue)],
  # num = num[lower.tri(cue)],
  # let = let[lower.tri(cue)],
  tar = tar[lower.tri(cue)],
  dis = dis[lower.tri(cue)],
  resp = resp[lower.tri(cue)],
  # iicc = incon[lower.tri(cue)]
  incon = incon[lower.tri(cue)]
  # congr = congr[is.lower.tri]
)
X <- 1 - X  ## convert to dissim

vifs <- car::vif(lm(1:nrow(X) ~ ., as.data.frame(X)))


dimnames(cue) <- list(regs, regs)
dimnames(tar) <- list(regs, regs)
dimnames(dis) <- list(regs, regs)
dimnames(resp) <- list(regs, regs)
dimnames(incon) <- list(regs, regs)


ctsmods <- list(cue = 1 - cue, tar = 1 - tar, dis = 1 - dis, resp = 1 - resp, incon = 1 - incon)  ## convert to dissim

## ----

saveRDS(ctsmods, here::here("out", "cuedts", "rsmods_full.RDS"))
data.table::fwrite(data.table::as.data.table(X), here::here("out", "cuedts", "rsmods_lt.csv"))