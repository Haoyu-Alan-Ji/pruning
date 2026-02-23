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
