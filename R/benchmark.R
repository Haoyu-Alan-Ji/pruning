library(corHMM)
library(rbenchmark)

## setname <- function(x) {cbind(nm = rownames(x), x)}
##' @param nstate number of states per trait (currently limited to 2)
##' @param ntrait number of traits
##' @param ntaxa number of taxa/phylogenty tips
##' @param seed random-number seed
##' @param meanrate mean transition rate
##' @param seql trait
##' @param collapse_allow for corHMM collapse function
##' @examples
##' simfun(ntaxa = 8)
simfun <- function(nstate = 2, ntrait = 2, ntaxa = 20, seed = NULL,
                   meanrate = 1, seql = 1, collapse = FALSE) {
  if (nstate!=2) stop("oops, simSeq to trait matrix not implemented for nstate!=2, 
                      to_bin's binary algorithm determines that the state can only be 2")
  if (!is.null(seed)) set.seed(seed)
  require("ape")
  require("phangorn")
  phy <- ape::rtree(ntaxa)
  phy <- reorder(phy, "pruningwise")
  
  Q <- Q_template(n=nstate, k= ntrait)
  nrates <- sum(Q != 0)
  Q[Q!=0] <- rexp(nrates, rate = 1/meanrate)
  
  s <- phangorn::simSeq(phy, l = seql, Q = Q,
                        type = "USER", levels = seq(nstate^ntrait),
                        rate = 1)
  if (collapse == FALSE) {
    repeat {
      s <- phangorn::simSeq(phy, l = seql, Q = Q,
                            type = "USER", levels = seq(nstate^ntrait),
                            rate = 1)
      if (nrow(unique(as.character(s))) == prod(nstate^ntrait)) break
    }
  }
  traitMatrix <- sapply(s, function(x) to_bin(x = x-1, n = ntrait)) |>
    t() |>
    as.data.frame() #|>
  #setname()
  
  list(tree = phy, data = traitMatrix)
}

parfun <- function(x) log(na.omit(c(x$solution)))

fitfun  <- function(dat, ..., rate.cat = 1) {
  tt <- system.time(
    invisible(capture.output(suppressWarnings(x <- with(dat, corHMM(tree, data, rate.cat = rate.cat, ...)))))
  )
  attr(x, "time") <- tt
  x
}

sumfun <- function(ntrait = 2, ntaxa = 200, model = "ARD", seed = NULL,
                    traitMatrix = NULL, realtree = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  seed <- seed %||% NA
  if (!is.null(traitMatrix) && !is.null(realtree)) {
  ss <- realfun(raw_traitM = traitMatrix,
                raw_tree = realtree) 
  } else {
  ss <- simfun(ntrait = ntrait, ntaxa = ntaxa, seed = seed)
  }

  fit_orig <- fitfun(ss, model = model, ...)
  fit_RTMB <- fitfun(ss, use_RTMB = TRUE, model = model, ...)
  
  p_orig <- parfun(fit_orig)
  p_RTMB <- parfun(fit_RTMB)
  p_diff <- p_orig-p_RTMB
  ## truncate par values at -10 (these diffs are mostly irrelevant)
  p_orig_trunc <- pmax(p_orig, -10)
  p_RTMB_trunc <- pmax(p_RTMB, -10)
  p_diff_trunc <- p_orig_trunc - p_RTMB_trunc
  ## get RMSE vs true rates on truncated scale ...
  if (!is.null(ss$true_rates)) {
    p_RTMB_rmse <- sqrt(mean((p_RTMB_trunc - ss$true_rates)^2))
    p_orig_rmse <- sqrt(mean((p_orig_trunc - ss$true_rates)^2))
  } else {
    p_RTMB_rmse <- NA
    p_orig_rmse <- NA
  }
  data.frame(seed, ntrait, ntaxa, model,
             RTMB_opt.time = fit_RTMB$opt.time[["elapsed"]],
             orig_opt.time = fit_orig$opt.time[["elapsed"]],
             RTMB_tot.time = attr(fit_RTMB, "time")[["elapsed"]],
             orig_tot.time = attr(fit_orig, "time")[["elapsed"]],
             RTMB_loglik = fit_RTMB$loglik,
             orig_loglik = fit_orig$loglik,
             par.rmsdiff = sqrt(mean(p_diff^2)),
             par.maxdiff = max(abs(p_diff)),
             par.rmsdiff.trunc = sqrt(mean(p_diff_trunc^2)),
             par.maxdiff.trunc = max(abs(p_diff_trunc)),
             RTMB_rmse = p_RTMB_rmse,
             orig_rmse = p_orig_rmse
  )
}

## convenience function to get everything we need for examples
## FIXME: make everything into an R package, this would all happen automatically ...
get_all <- function() {
    eval.parent({
        sourcefiles <- c("postAD.R", "general.R", "realtree.R", "benchmark.R")
        invisible(sapply(here::here(sprintf("R/%s", sourcefiles)), source))
        assign("c_traits", readRDS(here::here("data/c_traits.rds")), parent.frame())
        assign("c_tree", readRDS(here::here("data/c_tree.rds")), parent.frame())
    })
    invisible(NULL)
}
    
##' @examples
##' get_all()
##' dat <- list(tree = c_tree, data = c_traits)
##' system.time(pp0 <- adfit(dat = dat, keep_all = TRUE))
##' flist <- list(ag~care*spawning, care~1, spawning~1)
##' system.time(pp1 <- adfit(dat = dat, formula_list = flist, keep_all = TRUE))
##' stopifnot(length(pp0$pars.best) == 24)
##' stopifnot(length(pp1$pars.best) == 12)
##' nll_diff <- pp0$obj.best - pp1$obj.best ## NLL difference, rates minus formula
##' ## likelihood ratio test
##' pchisq(-2*nll_diff, df = 12, lower.tail = FALSE)
##' ## trying to diagnose why NLL for unconstrained (rates) is *worse* (larger)
##' ##  than constrained (formula) when we use random data (doesn't make sense,
##' ##  because constrained is nested in unconstrained, so NLL for unconstrained
##' ##  should always be *lower* (better), although not much better, especially
##' ##  for random trait data! Can this be a getting-stuck/multi-start problem??
##' set.seed(101)
##' rdat <- realfun(2, c_traits, c_tree)
##' ppr0 <- adfit(dat = rdat, keep_all = TRUE, multistart = 10, parallel = TRUE)
##' ppr1 <- adfit(dat = rdat, formula_list = flist, keep_all = TRUE, multistart = 10,
##'     parallel = TRUE)
##' ppr0$obj.best - ppr1$obj.best
##' rf0 <- ppr0$result_frame
##' rf1 <- rf0[ , grepl("^fitted", colnames(rf0)) ]
##' pairs(rf1[,1:12], gap = 0)
adfit <- function(dat, ...) {
  tt <- system.time(invisible(capture.output(suppressWarnings(x <- with(dat, postAD(tree = tree, traitMatrix = data,...))))))
  attr(x, "time") <- tt
  x
}

##' @examples
##' get_all()
##' t2 <- adsum(seed = 101, traitMatrix = c_traits, realtree = c_tree, formula_list = list(ag~care*spawning, care~1, spawning~1))
adsum <- function(ntaxa = 200, state = 2, seed, traitMatrix, realtree, formula_list, ...) {
  if (!is.null(seed)) set.seed(seed)
  seed <- seed %||% NA
  ntrait <- ncol(traitMatrix) - 1L
  ss <- realfun(raw_traitM = traitMatrix, raw_tree = realtree) 
  
  fit_rates <- adfit(ss, state = state, keep_all = TRUE, ...)
  fit_formula <- adfit(ss, state = state, formula_list = formula_list, keep_all = TRUE, ...)

  data.frame(
    seed, ntrait, state = state, ntaxa,
    rates_opt.time = fit_rates$result_frame$time_opt,
    formula_opt.time = fit_formula$result_frame$time_opt,
    rates_tot.time = attr(fit_rates, "time")[["elapsed"]],
    formula_tot.time = attr(fit_formula, "time")[["elapsed"]],
    rates_loglik = fit_rates$obj.best,
    formula_loglik = fit_formula$obj.best,
    rates_grad_norm = sqrt(sum(fit_rates$gr.best^2)),
    formula_grad_norm = sqrt(sum(fit_formula$gr.best^2))
  )
}
