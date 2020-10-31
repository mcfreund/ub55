tidy_lmer <- function(x) {
  coefs <- as.data.frame(coef(summary(x)))
  coefs$term <- rownames(coefs)
  dplyr::rename(coefs, estimate = "Estimate", se = "Std. Error", "statistic" = "t value", p.value = "Pr(>|t|)")
}


read_betas <- function(
  .subjs,
  .task,
  .glm,
  .dir
) {
  # .subjs = subjs; .task = "Cuedts"; .glm = "baseline_CongruencySwitch_EVENTS_censored_shifted"; .dir = dir.analysis
  
  ## initialize array
  
  pick.a.file <- 
    file.path(.dir, .subjs[1], "RESULTS",  .task, paste0(.glm, "_1"), paste0("STATS_", subjs[1], "_1_L_REML.func.gii"))
  labs <- afni("3dinfo", paste0("-label ", pick.a.file))
  labs <- unlist(strsplit(labs, "\\|"))
  is.reg <- !grepl("Full|block|Tstat|Fstat", labs)
  tab <- do.call(rbind, strsplit(gsub("_Coef", "", labs[is.reg]), "#"))
  trs <- as.numeric(unique(tab[, 2])) + 1
  regs <- unique(tab[, 1])
  
  
  n.vertex <- 10242
  n.tr <- length(trs)
  n.reg <- length(regs)
  n.subj <- length(subjs)
  
  betas <- array(
    NA,
    dim = c(n.vertex*2, n.reg, n.tr, n.subj, 2),
    dimnames = list(vertex = NULL, reg = regs, tr = NULL, subj = subjs, run = c("run1", "run2"))
  )
  
  vertex.inds <- cbind(L = 1:n.vertex, R = (n.vertex + 1):(n.vertex * 2))
  
  for (subj.i in seq_along(subjs)) {
    # subj.i = 1; run.i = 1; hemi.i = "L"
    
    for (run.i in 1:2) {
      # run.i = 1
      
      for (hemi.i in c("L", "R")) {
        # hemi.i = "L"
        
        inds <- vertex.inds[, hemi.i]
        
        fname <- file.path(
          .dir, .subjs[subj.i], "RESULTS",  .task, paste0(.glm, "_", run.i),  
          paste0("STATS_", subjs[subj.i], "_", run.i, "_", hemi.i, "_REML.func.gii")
        )
        
        if (!file.exists(fname)) next
        
        B <- mikeutils::read_gifti2matrix(fname)[is.reg, ]
        
        is.ok.i <- isTRUE(all.equal(dim(B), c(n.reg * n.tr, n.vertex)))
        if (!is.ok.i) stop("mismatched beta array")
        
        
        for (reg.i in seq_len(n.reg)) {
          # reg.i = 1
          
          is.reg.i <- grepl(paste0("^", regs[reg.i]), labs[is.reg])
          B.reg.i <- t(B[is.reg.i, ])
          
          is.ok.ii <- isTRUE(all.equal(dim(betas[inds, reg.i, , subj.i, run.i]), dim(B.reg.i)))
          if (!is.ok.ii) stop("mismatched regressor array")
          
          betas[inds, reg.i, , subj.i, run.i] <- B.reg.i
          
          }
        
      }

    }
      
  }
  
  betas
    
}




symmat4ggplot <- function(R, var.names = c("v1", "v2"), val.name = "value") {
  
  ## make factors for row and column labels
  dn <- dimnames(R)
  if (is.null(dn)) {
    dn <- setNames(list(paste0("cell_", 1:nrow(R)), paste0("cell_", 1:ncol(R))), var.names)
  } else {
    names(dn) <- var.names  
  }
  
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  labels[[2]] <- factor(labels[[2]], levels = rev(levels(labels[[2]])))
  
  r <- c(R)
  
  cbind(labels, setNames(as.data.frame(c(R)), val.name))
  
}




read_resid <- function(
  .subj,
  .task,
  .glm,
  .run,
  .dir
) {
  # .subj = subjs[1]; .task = "Axcpt"; .glm = "baseline_Cues_EVENTS_censored_shifted"; .dir = dir.analysis; run.i = 1
  
  
  
  fname_l <- file.path(
    .dir, .subj, "RESULTS",  .task, 
    paste0(.glm, "_", .run),
    paste0("wherr_", .subj, "_", .run, "_L_REML.func.gii")
  )
  fname_r <- gsub("_L_REML", "_R_REML", fname_l)
  
  E_l <- mikeutils::read_gifti2matrix(fname_l)
  E_r <- mikeutils::read_gifti2matrix(fname_r)
  
  cbind(E_l, E_r)
  
}


get.network <- function(x) {
  gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
}


loads <- function(x, dims = c("PC1", "PC2")) {
  
  rot <- as.data.frame(x$rotation)[, dims]
  rot$variable <- rownames(rot)
  
  rot
  
}


distance_cv <- function(B1, B2, m, regressors) {
  # B1 = B[, , 1]; B2 = B[, , 2]
  
  D <- colSums(t(m %*% B1 * m %*% B2)) / ncol(B1)
  matrix(D, ncol = length(regressors), dimnames = list(.row = regressors, .col = regressors))
  
}

enlist <- function(x) setNames(vector("list", length(x)), x)



subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')
  
  sub_chunk <- paste0("
  `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(", 
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}




matplot <- function(x) {
  
  ggplot(symmat4ggplot(x), aes(v1, v2, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "inferno") +
  theme_minimal() +
  theme(
    axis.text = element_blank(), axis.title = element_blank(), legend.position = "none",
    panel.border = element_blank(), panel.grid = element_blank()
    )
  
}



pdist2 <- function(A,B) {
  
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  
  tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
}



fprint <- function(B1, B2) {
  
  ## dimensions: rows == subjects, colums == features
  
  if (!identical(dim(B1), dim(B2))) stop('mats have same dims')
  
  m <- ncol(B1)
  n <- nrow(B1)
  
  D <- pdist2(B1, B2) / m  ## pairwise squared euclidean distances, per feature
  
  d_bn <- (colSums(D) + rowSums(D) - 2*diag(D)) / ((n-1)*2)  ## mean between-subj distances
  ##...colSums corresponds to using run2 as "template" and run1 as "database"
  ##...rowSums corresponds to using run1 as "template" and run2 as "database"
  ##...however 2*diagonal needs to be removed from those sums (as diagonal included in both colsums and rowsums.)
  ##...then we need to scale by twice the number of btw-subj comparisons per subject, (nsubjs-1)*2, to get the mean.
  
  contrast <- d_bn - diag(D)  ## bn-subj dists minus within subject dists.
  ##...positive contrast indicates consistent, subject-specific idiosyncracy in effect of condition.
  
  list(contrast = contrast, D = D)
  
}






## spot-checking functions:
# run1subj <- 25
# run2subj <- 24
# 
# 
# (pdistfun <- (res.mv$D[run1subj, run2subj] * nrow(B1))^2)
# (manual <- sum((B1[run1subj, ] - B2[run2subj, ])^2))
# (distfun <- c(dist(rbind(B1[run1subj, ], B2[run2subj, ]), method = "euclidean")^2))
# all.equal(manual, distfun)
# all.equal(manual, pdistfun)
# 
# (manual.uv <- (B1_bar[run1subj, ] - B2_bar[run2subj, ])^2)
# (pdistfun.uv <- (res.uv$D[run1subj, run2subj]))
# all.equal(manual.uv[[1]], pdistfun.uv)
# 
# all.equal(mean(c(D[1, 2:54], D[2:54, 1])), d_bn[[1]])
# D[1, 1]
# mean(c(D[1, 2:54], D[2:54, 1])) - D[1, 1]

