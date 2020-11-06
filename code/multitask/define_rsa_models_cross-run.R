
source(here::here("code", "_packages.R"))
source(here("code", "read-behav.R"))
source(here("code", "_vars.R"))
suppressWarnings(source(here("code", "_atlases.R")))
source(here("code", "_settings.R"))
source(here("code", "_funs.R"))


## base factors ----

conditions <- combo_paste(tasks, c("hi", "lo"), c("run1", "run2"))
n.dim <- length(conditions)

rsm.empty <- matrix(0, ncol = n.dim, nrow = n.dim, dimnames = list(conditions, conditions))
rsm.run <- rsm.empty
rsm.task <- rsm.empty
rsm.hi <- rsm.empty
rsm.lo <- rsm.empty

for (task.i in tasks) rsm.task[grepl(task.i, conditions), grepl(task.i, conditions)] <- 1
for (run.i in c("run1", "run2")) rsm.run[grepl(run.i, conditions), grepl(run.i, conditions)] <- 1
rsm.hi[grepl("hi", conditions), grepl("hi", conditions)] <- 1
rsm.lo[grepl("lo", conditions), grepl("lo", conditions)] <- 1
# diag(rsm.hi) <- 1
# diag(rsm.lo) <- 1


rsm.run %>% matplot + theme(axis.text.y = element_text()) + labs(title = "run")
rsm.task %>% matplot + theme(axis.text.y = element_text()) + labs(title = "task")
rsm.hi %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-high")
rsm.lo %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-low")
(1 - (rsm.run | rsm.task | rsm.hi | rsm.lo)) %>%
  matplot + theme(axis.text.y = element_text()) + labs(title = "intercept")




## within and between-run models ----

rsm.run.w <- rsm.run
rsm.run.b <- 1 - rsm.run
rsm.task.w <- rsm.run.w * rsm.task
rsm.task.b <- rsm.run.b * rsm.task
rsm.hi.w <- rsm.run.w * rsm.hi
rsm.hi.b <- rsm.run.b * rsm.hi
rsm.lo.w <- rsm.run.w * rsm.lo
rsm.lo.b <- rsm.run.b * rsm.lo
rsm.b0.w <- (1 - (rsm.task | rsm.hi | rsm.lo)) * rsm.run.w
rsm.b0.b <- (1 - (rsm.task | rsm.hi | rsm.lo)) * rsm.run.b


## check 


rsm.run.w  %>% matplot + theme(axis.text.y = element_text()) + labs(title = "within-run")
rsm.run.b  %>% matplot + theme(axis.text.y = element_text()) + labs(title = "between-run")

rsm.task.w %>% matplot + theme(axis.text.y = element_text()) + labs(title = "task within-run")
rsm.hi.w   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-high within-run")
rsm.lo.w   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-low within-run")
rsm.b0.w   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "intercept within-run")
rsm.task.b %>% matplot + theme(axis.text.y = element_text()) + labs(title = "task between-run")
rsm.hi.b   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-high between-run")
rsm.lo.b   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "condition-low between-run")
rsm.b0.b   %>% matplot + theme(axis.text.y = element_text()) + labs(title = "intercept between-run")


## melt to data.frame

tris <- cbind(
  
  
  mat2vec(rsm.task.w, value.name = "w_task"),
  w_hi = rsm.hi.w[lower.tri(rsm.empty)],
  w_lo = rsm.lo.w[lower.tri(rsm.empty)],
  w_b0 = rsm.b0.w[lower.tri(rsm.empty)],
  
  b_task = rsm.task.b[lower.tri(rsm.empty)],
  b_hi = rsm.hi.b[lower.tri(rsm.empty)],
  b_lo = rsm.lo.b[lower.tri(rsm.empty)],
  b_b0 = rsm.b0.b[lower.tri(rsm.empty)]
  
)


## examine xmat

X <- tris %>% select(where(is.numeric)) %>% as.matrix

cor(X)  ## good

kappa(X)  ## good
car::vif(lm(rep(1, nrow(tris)) ~ . + 0, data = as.data.frame(X)))  ## good

rowSums(X)  ## all modeled



## write ----

fwrite(tris, here("out", "multitask", "rsa-model-tris_cross-run.csv"))


