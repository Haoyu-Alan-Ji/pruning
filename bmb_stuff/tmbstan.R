setwd(here::here())
source("R/pruning_funs.R")
source("R/postAD.R")
library(tmbstan)

set.seed(427)
m3 <- rtree(20)
g1 <- reorder(m3, "pruningwise")
log_trans_rates <- log(abs(rnorm(8)))
objfun <- postAD(g1, 2, 2, log_trans_rates, multistart = 10, seed = 427, return_obj = TRUE)
objfun$fn(objfun$par)

options(mc.cores = 8)
stanfit <- tmbstan(objfun, chains=8, iter = 8000, seed = 20260302)

## https://bbolker.github.io/bbmisc/bayes/

## diagnostics: R-hat, ESS, etc
library(bayestestR)
diagnostic_posterior(stanfit)

diagnostic_posterior(stanfit1)
## trace plots, ...
library(bayesplot)
mcmc_trace(stanfit)
mcmc_rank_overlay(stanfit)

## don't like this one ...
## m <- mcmc_dens_chains(stanfit, regex_pars = "log_.*rate")

## (univariate) marginal densities ...
mcmc_dens_overlay(stanfit, regex_pars = "log_.*rate")

## don't really like this, pairs() method is better (see below)
## mcmc_pairs(stanfit, regex_pars = "log_.*rate")

## library(shinystan)
## launch_shinystan(stanfit)

## bivariate marginal densities ...
pairs(stanfit, gap = 0)  ## uses accept_prob to distinguish lower/upper triangle
try(pairs(stanfit, condition = "divergent__", gap = 0))
pairs(stanfit, gap = 0, condition = "stepsize__")

## now try it with ray-finned fish example!
## or ... could also try it (maybe first??) with some of the built-in corHMM
##  examples, e.g. the primate example ...
