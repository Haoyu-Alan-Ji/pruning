library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)
library(future)
library(future.apply)
library(parallel)

setwd(here::here())
source("R/ADtools.R")

pack_refit <- function(start, res) {
  start <- as.numeric(start)
  fitted <- as.numeric(res$par)

  names(start)  <- paste0("start",  seq_along(start))
  names(fitted) <- paste0("fitted", seq_along(fitted))

  data.frame(
    objective   = res$objective,
    convergence = res$convergence,
    message     = res$message,
    rbind(start),
    rbind(fitted),
    check.names = FALSE
  )
}

random_refit <- function(task, Phylodata, pars0, opt.args = NULL) {

  t0 <- proc.time()[["elapsed"]]

  ff <- RTMB::MakeADFun(
    func = cmb(postfun, Phylodata),
    parameters = list(log_trans_rates = pars0),
    silent = TRUE
  )

  t1 <- proc.time()[["elapsed"]]

  res <- with(ff, do.call(nlminb, c(list(task$start, fn, gr), opt.args)))

  t2 <- proc.time()[["elapsed"]]

  df <- pack_refit(task$start, res)

  df$refit_id    <- task$id
  df$start_seed  <- task$seed
  df$pid         <- Sys.getpid()
  df$nodename    <- Sys.info()[["nodename"]]
  df$time_MakeADFun   <- t1 - t0          # MakeADFun
  df$time_opt    <- t2 - t1          # nlminb 
  df$time_total  <- t2 - t0   

  df
}

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

postAD_parallel <- function(tree, trait, state, pars,
                            traitMatrix = NULL,
                            multistart = 10,
                            parallel = TRUE,
                            seed = 427,
                            jitter.sd = 0.5,
                            opt.args = NULL,
                            keep_all = FALSE,
                            rng_misuse = c("warning","error","ignore")) {

  rng_misuse <- match.arg(rng_misuse)

  t_total0 <- proc.time()[["elapsed"]]

  set.seed(seed)
  tasks <- lapply(seq_len(multistart - 1), function(i) {
    si <- seed + i
    set.seed(si)
    list(
      id   = i,
      seed = si,
      start = as.numeric(pars) + rnorm(length(pars), 0, jitter.sd)
    )
  })

  Q_temp <- Q_template(state, trait)

  if (!is.null(traitMatrix)) {
    repeat {
      set.seed(seed)
      s <- phangorn::simSeq(tree, l = 1, Q = Q_temp,
                            type = "USER", levels = seq(state^trait), rate = 1)
      if (nrow(unique(as.character(s))) == prod(state^trait)) break
    }
    s <- as.numeric(unlist(s))
  } else {
    multi_to_single <- function(traits, n = NULL) {
      x <- rev(cumprod(rev(n)))
      x <- c(x[-1], 1) 
      rowSums(sweep(traits, MARGIN = 2, x, "*")) + 1
    }
    s <- multi_to_single(traitMatrix, c(trait, state))
  }

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
  Phylodata <- list(Q_template = Q_temp, tree = tree, trait_values = s,
                    gainloss_pairs = gainloss_pairs)

  # baseline
  t5 <- proc.time()[["elapsed"]]
  ff0 <- RTMB::MakeADFun(
    func       = cmb(postfun, Phylodata),
    parameters = list(log_trans_rates = pars),
    silent     = TRUE
  )

  t6 <- proc.time()[["elapsed"]]
  res0 <- with(ff0, do.call(nlminb, c(list(ff0$par, ff0$fn, ff0$gr), opt.args)))
  t7 <- proc.time()[["elapsed"]]
  df0  <- pack_refit(ff0$par, res0)
  df0$refit_id   <- 0L
  df0$start_seed <- seed
  df0$pid        <- Sys.getpid()
  df0$nodename   <- Sys.info()[["nodename"]]
  df0$time_MakeADFun  <- t6 - t5
  df0$time_opt   <- t7 - t6 #nlminb
  df0$time_total <- t7 - t5
  
  workers <- 1
  if (parallel && multistart > 1) {
    workers <- min(multistart - 1, parallel::detectCores() - 1)
    workers <- max(workers, 1)

    # plan closure
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (workers > 1) {
      future::plan(future::cluster, workers = workers)
    } else {
      future::plan(future::sequential)
    }

    old_opt <- getOption("future.rng.onMisuse")
    on.exit(options(future.rng.onMisuse = old_opt), add = TRUE)
    options(future.rng.onMisuse = rng_misuse)
  }
    
  out_list <- if (multistart > 1) {
    future.apply::future_lapply(tasks, random_refit, Phylodata = Phylodata, pars0 = pars,
                                opt.args = opt.args, future.seed = TRUE
                                )
  } else {
    list()
  }
  result_frame <- do.call(rbind, c(list(df0), out_list))

  result_good <- subset(result_frame, convergence == 0)
  if (nrow(result_good) == 0) stop("No successful refits (convergence==0).")

  K <- length(ff0$par)
  fit_cols <- grep("^fitted", names(result_good), value = TRUE)[seq_len(K)]
  best_row <- which.min(result_good$objective)
  pars.best <- as.numeric(result_good[best_row, fit_cols])

  obj.best <- ff0$fn(pars.best)
  gr.best  <- ff0$gr(pars.best)

  t_total1 <- proc.time()[["elapsed"]]

  out <- list(pars.best  = pars.best, obj.best = obj.best, gr.best = gr.best, Phylo = Phylodata,
              multistart = multistart, workers = workers, cores_detected = parallel::detectCores(),
              time_total = t_total1 - t_total0
              )
  if (keep_all) out$result_frame <- result_frame
  out
}

res <- postAD_parallel(g1, 2, 2, log_trans_rates, multistart = 100, seed = 427)
res$obj.best
res$pars.best
