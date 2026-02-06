setwd(here::here())
library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)

set.seed(427)

source(here::here("R", "bestfit.R"))

## disabled opt.time and the output is a list
sum_test <- function(ntrait = 2, ntaxa = 200, model = "ARD", multistart = FALSE, 
                   try.times = 10, jitter.sd = 0.25, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  seed <- seed %||% NA
  ss <- simfun(ntrait = ntrait, ntaxa = ntaxa, seed = seed)
  if (multistart == TRUE) {
    fit_orig <- bestfit(ss, model = model, 
                        multistart = try.times, jitter.sd = jitter.sd, ...)
    fit_RTMB <- bestfit(ss, use_RTMB = TRUE, model = model, 
                        multistart = try.times, jitter.sd = jitter.sd, ...)
  } else {
    fit_orig <- fitfun(ss, model = model, ...)
    fit_RTMB <- fitfun(ss, use_RTMB = TRUE, model = model, ...)
    }
  p_orig <- parfun(fit_orig)
  p_RTMB <- parfun(fit_RTMB)
  p_diff <- p_orig-p_RTMB
  ## truncate par values at -10 (these diffs are mostly irrelevant)
  p_orig_trunc <- pmax(p_orig, -10)
  p_RTMB_trunc <- pmax(p_RTMB, -10)
  p_diff_trunc <- p_orig_trunc - p_RTMB_trunc
  ## get RMSE vs true rates on truncated scale ...
  p_RTMB_rmse <- sqrt(mean((p_RTMB_trunc - ss$true_rates)^2))
  p_orig_rmse <- sqrt(mean((p_orig_trunc - ss$true_rates)^2))
  list(seed, ntrait, ntaxa, model,
             ## RTMB_opt.time = fit_RTMB$opt.time[["elapsed"]],
             ## orig_opt.time = fit_orig$opt.time[["elapsed"]],
             RTMB_tot.time = attr(fit_RTMB, "time")[["elapsed"]],
             orig_tot.time = attr(fit_orig, "time")[["elapsed"]],
             RTMB_loglik = fit_RTMB$loglik,
             orig_loglik = fit_orig$loglik,
             par.rmsdiff = sqrt(mean(p_diff^2)),
             par.maxdiff = max(abs(p_diff)),
             par.rmsdiff.trunc = sqrt(mean(p_diff_trunc^2)),
             par.maxdiff.trunc = max(abs(p_diff_trunc)),
             RTMB_rmse = p_RTMB_rmse,
             orig_rmse = p_orig_rmse,
             fit_orig = fit_orig,
             tit_RTMB = fit_RTMB
  )
}

o2 <- sum_test(seed = 105, ntrait = 3, ntaxa = 20, multistart = TRUE, try.times = 20)
o1 <- sumfun(seed = 105, ntrait = 3, ntaxa = 20, multistart = TRUE, try.times = 20)

o2$fit_orig$opt.time
## NULL

o2$tit_RTMB$opt.time
## NULL

class(o1)
## [1] "data.frame"

dd <- readRDS(here::here('data', 'multi1.rds'))
class(dd)
## [1] "matrix" "array" 