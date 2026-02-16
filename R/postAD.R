library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)

setwd(here::here())

source("R/ADtools.R")

#' @param p parameters (log-hazard rates)
#' @param lb lower bound(s) for baseline priors
#' @param ub upper bound(s)
#' @param range width of Gaussian (+/- SD between mean and lower/upper bounds)
#' @param gainloss_pairs
#' @param lb_gainloss
#' @param ub_gainloss
#' @param range_gainloss number of SDs from center to lower/upper bounds
#' @param nllfun \emph{negative} log-likelihood function
#' @param negative return negative log posterior?


postfun <- function(pars, Phylodata,
                    ## add whatever arguments the RTMB pruning algorithm loglik function
                    ## needs (tree, trait data, etc.)
                    # p,
                    lb = log(1e-9), ub = log(1e2), range = 3,
                    # gainloss_pairs = NULL,
                    lb_gainloss = log(1e-3), ub_gainloss = log(1e3), range_gainloss = 3,
                    negative = FALSE
                    ) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)
  
  ## call RTMB pruning-algorithm code here to compute log-likelihood ...
  nll <- prune_nll(pars, Phylodata)
  loglik <- -1*nll

  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  p <- log_trans_rates # should I do this?
  logdnorm <- function(x, mu, sd) {
    -0.5 * log(2*pi*sd^2) - 0.5 * ((x - mu)/sd)^2
  }
  log.prior  <- sum(logdnorm(p, prior.mean, prior.sd))
  ## product of likelihood and prior -> sum of LL and log-prior
  res <- loglik + log.prior
  ## calculate gain/loss priors

  gl.prior.mean <- (lb_gainloss + ub_gainloss) / 2
  gl.prior.sd <- (ub_gainloss - lb_gainloss) / (2 * range_gainloss)
  ## vapply might not work in RTMB? replace with for loop?
  gl.values <- p[1] * 0 + numeric(length(gainloss_pairs))
  for (i in seq_len(length(gainloss_pairs))) {
    idx <- gainloss_pairs[[i]]
    gl.values[i] <- p[idx[1]] - p[idx[2]]
  }
  gl.log.prior <- sum(logdnorm(gl.values, gl.prior.mean, gl.prior.sd))

  res <- -1*(res + gl.log.prior)
  
  return(res)
}


postAD <- function(tree, trait, state, pars, rep.times = 100, opt.args = NULL) {
  
  Q_temp <- Q_template(state,trait)

  repeat {
  s <- phangorn::simSeq(tree, l = 1, Q = Q_temp,
                        type = "USER", levels = seq(state^trait),
                        rate = 1)
  if (nrow(unique(as.character(s))) == prod(state^trait)) break
  }
  s <- as.numeric(unlist(s))

  gl_pairs <- function(Qtemp) {
    if (!is.matrix(Qtemp) || nrow(Qtemp) != ncol(Qtemp)) {
      stop("Qtemp must be a square matrix.")
    }

    idx <- which(upper.tri(Qtemp) & Qtemp != 0, arr.ind = TRUE)

    pairs <- vector("list", nrow(idx))
    k <- 0L

    for (r in seq_len(nrow(idx))) {
      i <- idx[r, "row"]
      j <- idx[r, "col"]

      a <- Qtemp[i, j]
      b <- Qtemp[j, i]

      if (a != 0 && b != 0) {
        k <- k + 1L
        pairs[[k]] <- c(max(a, b), min(a, b)) 
      }
    }
    return(pairs)
  }
  gainloss_pairs <- gl_pairs(Q_temp)

  Phylodata <- list(Q_template = Q_temp, tree = tree, trait_values = s, gainloss_pairs = gainloss_pairs)
  
  ff <- MakeADFun(cmb(postfun, Phylodata), list(log_trans_rates = pars), silent =TRUE)
  AD <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr), opt.args))))
  
  result_frame <- replicate(rep.times, simplify = FALSE, random_refit(ff)) |>
    do.call(what = rbind)
  
  result_good <- subset(result_frame, convergence == 0)

  fit_cols <- grep("^fitted", names(result_good), value = TRUE)
  fit_cols <- fit_cols[seq_len(length(ff$par))] 

  pars.best <- as.numeric(result_good[which.min(result_good$objective), fit_cols])
  
  gr.best <- ff$gr(pars.best)
  obj.best <- ff$fn(pars.best)
  fn.cdf <- ecdf(result_good$objective)
  
  result <- list(pars.best = pars.best,
                 obj.best = obj.best,
                 gr.best = gr.best,
                 Phylo = Phylodata,
                 fn.cdf = fn.cdf)
  
  cat('loglik:\n'); print(obj.best)
  cat('pars:\n'); print(pars.best)
  cat('loglik cdf:\n'); plot(fn.cdf)
  return(result)
}

set.seed(427)
m3 <- rtree(20)
g1 <- reorder(m3, "pruningwise")

log_trans_rates <- log(abs(rnorm(sum(Q_template4 != 0))))
t2 <- postAD(g1, 2, 2, log_trans_rates)
# loglik:
# [1] 52.11188
# pars:
# [1] -4.3221751  0.3149425 -3.1975299 -8.0608872 -0.3005529 -3.9103304 -8.0598451 -2.5813955