knitr::opts_chunk$set(
  cache = TRUE, echo = FALSE, warning = FALSE, message = FALSE,
  fig.align = 'center',
  cache.lazy = FALSE,  ##https://stackoverflow.com/questions/39417003/long-vectors-not-supported-yet-error-in-rmd-but-not-in-r-script
  # fig.width = 11.5,
  fig.fullwidth = TRUE
)

set.seed(0)

# theme_set(theme_minimal(base_size = 12))
theme_set(theme_half_open())
