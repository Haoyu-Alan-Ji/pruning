library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)
library(future)
library(future.apply)
library(parallel)

source(here::here('R', 'general.R'))

#' translate trait matrix (no species name!) to single trait
#' @param traits trait matrix (trait values, 0-indexed)
#' @param n number of states per trait
multi_to_single <- function(traits, n = NULL) {
  if (is.null(n)) {
      n <- apply(traits, 2, max)+1
      warning("guessing number of traits per state from max()+1")
  }
  x <- rev(cumprod(rev(n)))
  x <- c(x[-1], 1) 
  rowSums(sweep(traits, MARGIN = 2, x, "*")) + 1
}

gl_pairs <- function(q_prep) {
  if (q_prep$mode == "rates") {
    Qind <- q_prep$Q_indicator
    idx <- which(upper.tri(Qind) & Qind != 0, arr.ind = TRUE)

    pairs <- vector("list", 0)
    for (r in seq_len(nrow(idx))) {
      i <- idx[r, "row"]
      j <- idx[r, "col"]

      a <- Qind[i, j]
      b <- Qind[j, i]

      if (a != 0 && b != 0) {
        pairs[[length(pairs) + 1L]] <- c(max(a, b), min(a, b))
      }
    }
    return(pairs)
  } else if (q_prep$mode == "formula") {
    nm <- q_prep$par_names
    gain_idx <- grep("_gain:", nm, fixed = TRUE)

    pairs <- vector("list", 0)
    for (i in gain_idx) {
      partner <- sub("_gain:", "_loss:", nm[i], fixed = TRUE)
      j <- match(partner, nm)
      if (!is.na(j)) {
        pairs[[length(pairs) + 1L]] <- c(i, j)
      }
    }
    return(pairs)
  }
}

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
  p <- q_par
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

pack_refit <- function(start, res) {
  start <- as.numeric(start)
  fitted <- as.numeric(res$par)

  names(start)  <- paste0("start",  seq_along(start))
  names(fitted) <- paste0("fitted", seq_along(fitted))

  data.frame(objective = res$objective, convergence = res$convergence, message = res$message,
             rbind(start),rbind(fitted), check.names = FALSE)
}

random_refit <- function(task, Phylodata, pars0, opt.args = NULL) {

  t0 <- proc.time()[["elapsed"]]

  ff <- RTMB::MakeADFun(
    func = cmb(postfun, Phylodata),
    parameters = list(q_par = pars0),
    silent = TRUE
  )

  t1 <- proc.time()[["elapsed"]]

  res <- with(ff, do.call(nlminb, c(list(task$start, fn, gr), opt.args)))

  t2 <- proc.time()[["elapsed"]]

  df <- pack_refit(task$start, res)

  df$refit_id <- task$id
  df$start_seed  <- task$seed
  df$pid <- Sys.getpid()
  df$nodename <- Sys.info()[["nodename"]]
  df$time_MakeADFun <- t1 - t0          # MakeADFun
  df$time_opt <- t2 - t1          # nlminb 
  df$time_total <- t2 - t0   

  df
}

#' will simulate random trait data if `traitMatrix` is missing
#' @param verbose print output?
#' @examples
#' set.seed(427)
#' g1 <- reorder(ape::rtree(20), "pruningwise")
#' log_trans_rates <- log(abs(rnorm(8)))
#' res <- postAD(tree = g1, trait = 2, state = 2, pars = log_trans_rates, multistart = 10, seed = 427)
#' objfun <- postAD(g1, 2, 2, log_trans_rates, multistart = 10, seed = 427, return_obj = TRUE)
#' objfun$fn(objfun$par)
## FIXME: trait, state should probably be 'ntrait', 'nstate'
postAD <- function(tree, trait, state, pars, formula_list = NULL, traitMatrix = NULL,
                  multistart = 10, parallel = TRUE, jitter.sd = 0.5,
                  seed = 427, rng_misuse = c("warning","error","ignore"),
                  opt.args = NULL, keep_all = FALSE,
                  return_obj = FALSE,
                  return_obj_type = c("AD", "raw"),
                  verbose = FALSE) {

  return_obj_type <- match.arg(return_obj_type)
  mode <- if (is.null(formula_list)) 'rates' else  'formula'
  rng_misuse <- match.arg(rng_misuse)
  t_total0 <- proc.time()[["elapsed"]]
  set.seed(seed)

  tasks <- lapply(seq_len(multistart - 1), function(i) {
    si <- seed + i
    set.seed(si)
    list(id   = i, seed = si, start = as.numeric(pars) + rnorm(length(pars), 0, jitter.sd))
  })

  q_prep <- Q_prep(mode = mode, formula_list = formula_list, nstate = state, ntrait = trait)

  if (is.null(traitMatrix)) {
    repeat {
      set.seed(seed)
      s <- phangorn::simSeq(tree, l = 1, Q = q_prep$Q_indicator,
                            type = "USER", levels = seq(state^trait), rate = 1)
      if (nrow(unique(as.character(s))) == prod(state^trait)) break
    }
    s <- as.numeric(unlist(s))
  } else {
    nvec <- if (length(state) == 1) rep(state, trait) else state
    spnames <- NULL

    if (is.data.frame(traitMatrix) && ncol(traitMatrix) > 0) {
        first_col <- traitMatrix[[1]]

        if (is.character(first_col) || is.factor(first_col)) {
            spnames <- as.character(first_col)
            traitMatrix <- traitMatrix[, -1, drop = FALSE]
        }
    }

    traitMatrix <- as.data.frame(traitMatrix)
    traitMatrix[] <- lapply(traitMatrix, as.numeric)

    s <- multi_to_single(traitMatrix, nvec)

    if (!is.null(spnames)) {
      stopifnot(length(spnames) == ape::Ntip(tree))
      stopifnot(all(spnames == tree$tip.label))
      names(s) <- spnames
    }
  }

  gainloss_pairs <- gl_pairs(q_prep) 
  Phylodata <- list(q_prep = q_prep, tree = tree, trait_values = s, gainloss_pairs = gainloss_pairs)

  if (return_obj  && return_obj_type == "raw") return(list(fn = postfun, pars = list(q_par = pars), nll = function(par = pars) prune_nll(list(q_par = par), Phylodata), data = Phylodata))
    
  # baseline
  t5 <- proc.time()[["elapsed"]]
  ff0 <- RTMB::MakeADFun(func = cmb(postfun, Phylodata),
                          parameters = list(q_par = pars),
                          silent = TRUE
                        )
  if (return_obj) return(ff0)
  
  t6 <- proc.time()[["elapsed"]]
  res0 <- with(ff0, do.call(nlminb, c(list(ff0$par, ff0$fn, ff0$gr), opt.args)))
  t7 <- proc.time()[["elapsed"]]
  df0  <- pack_refit(ff0$par, res0)
  df0$refit_id <- 0L
  df0$start_seed <- seed
  df0$pid <- Sys.getpid()
  df0$nodename <- Sys.info()[["nodename"]]
  df0$time_MakeADFun <- t6 - t5
  df0$time_opt <- t7 - t6 #nlminb
  df0$time_total <- t7 - t5
  
  workers <- 1
  if (parallel && multistart > 1) {
    workers <- min(multistart - 1, parallel::detectCores() - 1)
    workers <- max(workers, 1)

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
  } else list()
  
  result_frame <- do.call(rbind, c(list(df0), out_list))

  result_good <- subset(result_frame, convergence == 0)
  if (nrow(result_good) == 0) stop("No successful refits (convergence==0).")

  fit_cols <- grep("^fitted", names(result_good), value = TRUE)[seq_len(length(ff0$par))]
  best_row <- which.min(result_good$objective)
  pars.best <- as.numeric(result_good[best_row, fit_cols])
  obj.best <- ff0$fn(pars.best)
  gr.best  <- ff0$gr(pars.best)
  fn.cdf <- ecdf(result_good$objective)

  t_total1 <- proc.time()[["elapsed"]]
  time_total = t_total1 - t_total0

  out <- list(pars.best  = pars.best, obj.best = obj.best, gr.best = gr.best, Phylo = Phylodata,
              multistart = multistart, workers = workers, cores_detected = parallel::detectCores(),
              time_total = time_total
              )
  if (keep_all) out$result_frame <- result_frame
  if (verbose) {
     cat('loglik:\n'); print(obj.best)
     cat('pars:\n'); print(pars.best)
     cat('time:\n'); print(time_total)
     cat('loglik cdf:\n'); plot(fn.cdf)
  }
  out
}

if (FALSE) {
  set.seed(101)
  library(corHMM)
  data(primates)
  prim_phy <- reorder(multi2di(primates[[1]]), "pruningwise")
  ## now reorder to match tip labels and drop names column
  prim_data <- primates[[2]][match(primates[[2]][,1], prim_phy$tip), -1]
  pp <- postAD(prim_phy, trait = 2, state = 2, pars = rep(-2, 6),
          formula_list = list(T1~1, T2 ~ T1))
}                       
## testing
if (FALSE) {
    set.seed(101)
    tree1 <- reorder(ape::rtree(40), "pruningwise")

    qprep1 <- Q_prep(mode = "rates", nstate = 2, ntrait = 3)
    pars1 <- setNames(rep(log(0.2), qprep1$n_par), qprep1$par_names)

    res <- postAD(tree1, 2, 2, pars1, multistart = 10, seed = 427)

    library(corHMM)
    data(primates)
    source("R/postAD.R")
    prim_phy <- reorder(multi2di(primates[[1]]), "pruningwise")
    ## now reorder to match tip labels and drop names column
    prim_data <- primates[[2]][match(primates[[2]][,1], prim_phy$tip), -1]
    pp <- postAD(prim_phy, trait = 2, state = 2, pars = rep(-2, 6),
           formula_list = list(T1~1, T2 ~ T1),
           return_obj = TRUE,
           return_obj_type = "raw")

    ## debug(pp$fn)

    with(pp, fn(pars, data))
         
}
