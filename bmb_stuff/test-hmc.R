library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(RTMB)
library(microbenchmark)
library(corHMM)
library(tmbstan)

source("../R/pruning_funs.R")
load("../data/primates.rda")
m4 <- rtree(100)

g2 <- reorder(m4, "pruningwise")
Q_template8 <- setup_Q_template(2,3)
log_trans_rates <- log(abs(rnorm(24)))
                       
d <- get_multitrait(c(2 ,2, 2), g2)
traits8 <- multi_to_single(d[,-1], c(2, 2, 2))

TMBdata <- list(Q_template = Q_template8,
                       tree = g2, trait_values = traits8)
ff4 <- RTMB::MakeADFun(func = prune_nll,
                      parameters = list(log_trans_rates = log_trans_rates), silent = TRUE)

bounds <- rep(25, length(ff4$par))
fit <- with(ff4, nlminb(par, fn, gr, lower = -1*bounds, upper = bounds))

summary(fit$par)
tt <- tmbstan(ff4, lower = -1*bounds, upper = bounds)

