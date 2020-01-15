#---
title: "Improving measures of pollinator performance: Pollen transfer quantity and quality determine plant reproduction"
author: "Jamie R. Stavert"
date: "15 January 2020"
#---

####################################
#load packages----
####################################

library(survival)
library(dplyr)
library(glmmTMB)
library(emmeans)

####################################
#read in data----
####################################

#pollen deposition and pollen tube growth data
survival <- read.csv(file = "data/survival-data.csv", header=TRUE, row.names=NULL)

#fruit set data
fruit.set <- read.csv(file = "data/fruit-set.csv", header=TRUE, row.names=NULL)

#seed set data
seed.set <- read.csv(file = "data/seed-count.csv", header=TRUE, row.names=NULL)

####################################
#survival models----
####################################

#categorical cox model
surv.fit <- coxph(Surv(time, valbin) ~ Treatment, data=survival)
summary(surv.fit)

#remove other treatment levels and run continous model on visits only
survival.cont <- survival[! survival$Treatment %in% c("HP","O","C"), ] %>% droplevels()
survival.cont$No.visits <- as.numeric(as.character(survival.cont$No.visits))

#fit the continuous cox model
surv.fit.cont <- coxph(Surv(time, valbin) ~ No.visits, data=survival.cont)#fit cox model
summary(surv.fit.cont)

####################################
#pollination success models----
####################################

#pollination success continuous model
poll.succ.cont <- glmmTMB(value ~ No.visits*variable + (1|uID/stigmanum),
                  data = survival.cont,
                  ziformula=~1,
                  family="nbinom1")

summary(poll.succ.cont)

#run emtrends to test for differences between slopes and if slopes are different from zero
em <- emtrends(poll.succ.cont, pairwise ~ variable, var="No.visits", adjust="FDR")
summary(em, infer=c(TRUE,TRUE),null=0, adjust="FDR")

#pollination success categorical model
poll.succ.cat <- glmmTMB(value ~ tc*variable + (1|uID/stigmanum),
                 data = survival,
                 ziformula=~1,
                 family="nbinom1")

summary(poll.succ.cat)

#run lsmeans pairwise comparisons
em.comps <- emmeans(poll.succ.cat, pairwise ~ tc|variable, adjust = "fdr", type = "response")
CLD(em.comps, Letters = letters, level = .95, adjust = "fdr")

####################################
#fruit set models----
#####################################

#run categorical model
fs.mod <- glmmTMB(final ~ Treatment + (1|year/farm/block),
               family = binomial,
               data=fruit.set)

summary(fs.mod)

#calculate pairwise comparisons
fs.mod <- emmeans(fs.mod, ~ Treatment, type = "response", level = .95, adjust="fdr")
CLD(fs.mod, Letters = letters, adjust="fdr")

#run continuous model
#remove non-visit treatments
fruit.set.cont <- fruit.set[! fruit.set$Treatment %in% c("HP","O","C"), ] %>% droplevels()
fruit.set.cont$No.visits <- as.numeric(as.character(fruit.set.cont$No.visits))
fruit.set.cont <- fruit.set.cont[fruit.set.cont$No.visits<11,]%>%droplevels()

#model and use emmeans for pairwise comparision between treatments for initial fruit set
fs.cont.mod <- glmmTMB(final ~ No.visits + (1|year/farm/block),
                       data = fruit.set.cont,
                       family="binomial")
summary(fs.cont.mod)

####################################
#seed set models----
####################################

#run categorical model
ss.mod <- glmmTMB(seeds ~ treatment + (1|year/farm/block),
              data = seed.set,
              ziformula=~1,
              family="nbinom2")

summary(ss.mod)

#run emmeans pairwise comparisons
comps.seed <- emmeans(ss.mod, pairwise ~ treatment, adjust = "fdr", type = "response")#include model, then the pairwise comparison (i.e. the interaction)
CLD(comps.seed, Letters = letters, level = .95, adjust = "fdr")

#run continuous model
#remove non-visit treatments
seed.set.cont <- seed.set[! seed.set$treatment %in% c("hp","o","c"), ] %>% droplevels()
seed.set.cont$no.visits <- as.numeric(as.character(seed.set.cont$no.visits))
seed.set.cont <- seed.set.cont[seed.set.cont$no.visits<11,] %>% droplevels()

ss.cont.mod = glmmTMB(seeds ~ no.visits + (1|year/farm/block),
                      data = seed.set.cont,
                      ziformula=~1,
                      family="nbinom2")

summary(ss.cont.mod)

####################################
#END
####################################