tidy_lmer <- function(x) {
  coefs <- as.data.frame(coef(summary(x)))
  coefs$term <- rownames(coefs)
  dplyr::rename(coefs, estimate = "Estimate", se = "Std. Error", "statistic" = "t value", p.value = "Pr(>|t|)")
}
