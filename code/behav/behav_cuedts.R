library(here)
library(knitr)
library(data.table)
library(dplyr)
library(gridExtra)
library(mikeutils)
library(ggplot2)
library(lme4)
library(nlme)

source(here("code", "read-behav.R"))

logit2prob <- function(x) exp(x) / (1 + exp(x))

cuedts <- as.data.table(cuedts)

## look ----

plot(cuedts$rt)

cuedts.rt <- cuedts[rt > 0 & acc == 1 & switch != ""]

grid.arrange(
  
  cuedts.rt %>%

    ggplot(aes(rt)) +
    geom_density(fill = "slateblue", alpha = 0.3) +
    
    labs(title = "all correct trials") +
    theme(legend.position = "none"),
  
  cuedts.rt %>%
    
    .[switch != ""] %>%
    
    ggplot(aes(rt)) +
    geom_density(aes(fill = interaction(switch, trial.type)), alpha = 0.3) +
    
    labs(title = "sequence*trialtype") +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme(legend.position = c(0.75, 0.5)),
  
  ncol = 2
  
)


cuedts.rt %>%
  
  ggplot(aes(sample = rt)) +
  stat_qq(alpha = 0.8, size = 1) +
  stat_qq_line(size = 0.25) +
  
  facet_wrap(vars(subj)) +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank()) +
  labs(
    title    = "QQ rt versus normal"
  )

cuedts.rt %>%
  
  ggplot(aes(x = rt, group = subj)) +
  geom_density(color = rgb(0, 0, 0, 0.5), size = 2)

## prelim model ----

## prelim search of model space

cuedts.rt$cresp <- as.factor(cuedts.rt$cresp)

m0 <- lmer(rt ~ trial.type + switch + (1 | subj), cuedts.rt)
summary(m0)
mint <- update(m0, . ~ . + trial.type:switch)
summary(mint)

m1 <- update(m0, . ~ . + cue + , control = lmerControl(optimizer = "bobyqa"))
summary(m1)

m2 <- update(m0, . ~ . + cue + cresp, control = lmerControl(optimizer = "bobyqa"))
summary(m2)

m3 <- lmer(rt ~ trial.type + switch + cue + cresp + (cresp | subj), cuedts.rt)
summary(m3)

m4 <- lmer(
  rt ~ trial.type + switch + cue + cresp + (cresp + cue | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
  )
summary(m4)

m5 <- lmer(
  rt ~ trial.type + switch + cue + (cue | subj) + (1 | cresp),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m5)

m6 <- lmer(
  rt ~ trial.type + switch + (trial.type | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m6)

m7 <- lmer(
  rt ~ trial.type + switch + (switch | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m7)


m8 <- lmer(
  rt ~ trial.type + switch + cue + cresp + (cue + cresp | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m8)

plot(ranef(m8)$subj)

m9 <- lmer(
  rt ~ trial.type + switch + cue + cresp + (1 | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
anova(m8, m9)

m10 <- lmer(
  rt ~ cresp * trial.type + cue + switch + (cue + cresp | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m10)

plot(m10)

m11 <- lmer(
  rt ~ cresp * trial.type * cue + switch + (cue + cresp | subj),
  cuedts.rt, 
  control = lmerControl(optimizer = "bobyqa")
)
summary(m11)


cuedts.rt$resid <- resid(m10)
cuedts.rt %>%
  
  ggplot(aes(sample = resid)) +
  stat_qq(alpha = 0.8, size = 1) +
  stat_qq_line(size = 0.25) +
  
  facet_wrap(vars(subj)) +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank()) +
  labs(
    title    = "QQ resid rt versus normal"
  )


## cresp2*trial.typei??

table(cuedts.rt[trial.type == "InCon", c("cue", "stimuli", "cresp")])
table(cuedts.rt[trial.type == "Con", c("cue", "stimuli", "cresp")])

## index finger more potent competitor?

cuedts.rt %>%
  
  group_by(cresp, trial.type, subj) %>%
  summarize(rt = mean(rt)) %>%
  tidyr::pivot_wider(names_from = c("trial.type", "cresp"), values_from = rt) %>%
  mutate(
    tt1 = InCon_1 - Con_1,
    tt2 = InCon_2 - Con_2,
    int = tt2 - tt1,
    me  = (tt2 + tt1)/2
    ) %>%
  select(-c(Con_1, InCon_1, Con_2, InCon_2)) %>%
  tidyr::pivot_longer(cols = c("tt1", "tt2", "int", "me"), values_to = "rt") %>%
  
  ggplot() +
  geom_density(aes(rt, color = name), size = 2) +
  scale_color_viridis_d() +
  theme(legend.position = c(0.75, 0.75))

cuedts.rt %>%
  
  group_by(cresp, trial.type, subj) %>%
  summarize(rt = mean(rt)) %>%
  
  ggplot(aes(interaction(cresp, trial.type), rt)) +
  stat_summary(fun.data = mean_cl_boot)


cuedts.rt %>%
  
  group_by(cresp, trial.type, subj, cue) %>%
  summarize(rt = mean(rt)) %>%

  ggplot(aes(interaction(cresp, cue, trial.type), rt)) +
  stat_summary(fun.data = mean_cl_boot)




## errors

cuedts$er <- 1 - cuedts$acc
cuedts <- cuedts[switch != ""]
cuedts$cresp <- as.factor(cuedts$cresp)

g0 <- glmer(
   er ~ trial.type + switch + (1 | subj),
  cuedts,
  family = binomial
  # control = lmerControl(optimizer = "bobyqa")
)
summary(g0)

gint <- update(g0, . ~ . + trial.type:switch)
summary(gint)

g1 <- update(g0, . ~ . + cue)
summary(g1)

g2 <- update(g0, . ~ . + cue + cresp)
summary(g2)

g3 <- glmer(
  er ~ trial.type + switch + (trial.type | subj),
  cuedts,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(g3)

g4 <- glmer(
  er ~ trial.type + switch + (switch | subj),
  cuedts,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(g4)


g5 <- glmer(
  er ~ cresp + trial.type + cue + switch + (1 | subj),
  cuedts,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(g5)


g6 <- glmer(
  er ~ cresp * trial.type + cue + switch + (1 | subj),
  cuedts,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(g6)
