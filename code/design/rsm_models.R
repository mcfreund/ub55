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



matplot <- function(x) {
  
  ggplot(symmat4ggplot(x), aes(v1, v2, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(option = "inferno") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), 
      legend.position = "none"
    )
  
}


model_matrix <- function(x, cols, l = NULL) {
  
  X <- model.matrix(~as.factor(x) - 1)
  colnames(X) <- gsub("as\\.factor\\(x\\)", "", colnames(X))
  
  XX <- tcrossprod(X[, cols])
  
  dimnames(XX) <- list(l, l)
  
  XX
  
}



library(ggplot2)
library(mikeutils)

conds <- combo_paste(c("letter", "number"), c("consonant", "vowel"), c("even", "odd"))
facts <- do.call(rbind, strsplit(conds, "_"))
colnames(facts) <- c("cue", "stimulus_l", "stimulus_n")

lt <- lower.tri(diag(length(conds)))


X_cue_letter <- model_matrix(facts[, "cue"], "letter", conds)
X_cue_number <- model_matrix(facts[, "cue"], "number", conds)
X_cue <- X_cue_number + X_cue_letter

X_stim_letter <- model_matrix(facts[, "stimulus_l"], c("consonant", "vowel"), conds)
X_stim_number <- model_matrix(facts[, "stimulus_n"], c("even", "odd"), conds)

X_target <- X_cue_letter*X_stim_letter + X_cue_number*X_stim_number
X_distractor <- X_cue_letter*X_stim_number + X_cue_number*X_stim_letter

X_vowel <- model_matrix(facts[, "stimulus_l"], "vowel", conds)
X_even <- model_matrix(facts[, "stimulus_n"], "even", conds)
X_consonant <- model_matrix(facts[, "stimulus_l"], "consonant", conds)
X_odd <- model_matrix(facts[, "stimulus_n"], "odd", conds)

X_incongruency <- tcrossprod(as.numeric(grepl("consonant_even|vowel_odd", conds)))
dimnames(X_incongruency) <- list(conds, conds)


matplot(X_cue[sort(conds), sort(conds)]) + labs(title = "cue")
matplot(X_target[sort(conds), sort(conds)]) + labs(title = "target")
matplot(X_distractor[sort(conds), sort(conds)]) + labs(title = "distractor")
matplot(X_incongruency[sort(conds), sort(conds)]) + labs(title = "incongruency")

v <- cbind(cue = X_cue[lt], target = X_target[lt], distractor = X_distractor[lt], incongruency = X_incongruency[lt])

matplot(cor(v)) + geom_text(aes(label = round(value, 2)), color = "grey40")

library("car")
car::vif(lm(rep(1, nrow(v)) ~ ., data.frame(v)))
