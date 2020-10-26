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
  # .subjs = subjs; .task = "Axcpt"; .glm = "baseline_Cues_EVENTS_censored_shifted"; .dir = dir.analysis
  
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
    # subj.i = 1
    
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
          
          is.reg.i <- grepl(regs[reg.i], labs[is.reg])
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
