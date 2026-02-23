get_multitrait <- function(trait_number, tree) {
  ss <- function(n) {
    sample(0:(n - 1), size = Ntip(tree), replace = TRUE)
  }
  repeat {
    d <- list()
    for (i in seq_along(trait_number)) {
      d[[i]] <- ss(trait_number[i])
    }
    d <- as.data.frame(d)
    names(d) <- paste0("trait", seq_along(trait_number))
    if (nrow(unique(d)) == prod(trait_number)) break
  }
  d <- cbind(tree$tip.label, d)
  return(d)
}

##' @examples
#' n <- c(2, 3, 5)
#' ## make up a random trait matrix
#' set.seed(101)
#' ss <- function(n) {
#'   sample(0:(n-1), size = 10, replace = TRUE)
#' }
#' d <- cbind(ss(2), ss(3), ss(5))
#' multitrait_to_int(d, c(2, 3, 5))
multi_to_single <- function(traits, n = NULL) {
  if (is.null(n)) {
    ## assume traits are zero-indexed
    n <- apply(traits, 2, max) + 1
  } else if (length(n) == 1) {
    n <- rep(n, ncol(as.data.frame(traits)))
  }
  x <- rev(cumprod(rev(n)))
  x <- c(x[-1], 1)  ## (15, 5, 1) for n = (2, 3, 5)
  rowSums(sweep(traits, MARGIN = 2, x, "*")) + 1
}

## TO CHECK: check against parameters? against full matrix?
#' @param corHMM_fit a fitted corHMM object
#' @examples
#' if (require("corHMM")) {
#' ## example from ?corHMM
#'     data(primates)
#'     phy <- multi2di(primates[[1]])
#'     data <- primates[[2]]
#'     MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
#'  ## try it in prune_nll
#'     TMBdata <- list(tree = phy,
#'          trait_values = multitrait_to_int(data[,-1], c(2,2)),
#'          Q_template = setup_Q_template(n=2, k=2))
#'     try(hmm2prune(MK_3state))
#' }
hmm2prune <- function(corHMM_fit) {
  index.mat <- corHMM_fit$index.mat
  solution <- corHMM_fit$solution
  if (ncol(TMBdata$Q_template) != ncol(index.mat)) {
    stop("incompatibility between data and corHMM matrix (trait collapse??)")
  }
  
  ids <- sort(na.omit(as.vector(index.mat)))
  trans_rates <- rep(0, length(ids))
  
  rate_list <- numeric(length(ids))
  for (i in seq_along(ids)) {
    id <- ids[i]
    p <- which(index.mat == id, arr.ind = TRUE)
    trans_rates[i] <- solution[p[1, 1], p[1, 2]]
  }
  pars <- list(log_trans_rates = log(trans_rates))
  result <- prune_nll(pars)
  return(result)
}

##' @example 
##' if (FALSE) {
##' system.time(Q_big <- setup_Q_template(k=5)) ## 10 seconds
##' png("bigmat.png")
##' imat(Q_big)
##' dev.off()
##' }
imat <- function(m, useRaster = TRUE, ...) {
  require(Matrix)
  image(Matrix(m), useRaster = useRaster, xlab = "", ylab = "", sub = "", ...)
}

## thought I could do a Kronecker product trick but now I don't see how ...
## if (use_kron) {
##   Q0 <- function(n) {
##     m <- matrix(1, n, n)
##     diag(m) <- 0
##     m
##   }
##   return(Reduce(kronecker, lapply(n, Q0)))
## }

# without prior
pruningAD <- function(tree, trait, state, rep.times = 100,
                      meanrate = 1, pars, opt.args = NULL) {
  Q <- Q_template(state,trait)

  repeat {
  s <- phangorn::simSeq(tree, l = 1, Q = Q,
                        type = "USER", levels = seq(nstate^ntrait),
                        rate = 1)
  if (nrow(unique(as.character(s))) == prod(nstate^ntrait)) break
  }

  s <- as.numeric(unlist(s))

  # nrates <- sum(Q != 0)
  # Q[Q!=0] <- rexp(nrates, rate = 1/meanrate)
  
  Phylodata <- list(Q_template = Q, tree = tree, trait_values = s)
  
  ff <- MakeADFun(cmb(prune_nll, Phylodata), list(log_trans_rates = pars), silent = TRUE)
  AD <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr), opt.args))))
  
  result_frame <- suppressWarnings(
    replicate(rep.times, simplify = FALSE, random_refit(ff))) |>
    do.call(what = rbind)
  
  result_good <- subset(result_frame, convergence == 0)

  if (nrow(result_good) == 0) {
    tt <- table(result_frame$convergence)
    stop("No good results found: all refits failed to converge: codes ",
         paste(sprintf("code %s = %d", names(tt), tt), collapse = ", "))
  }
  
  pars.best <- result_good[which.min(result_good$objective),
                          paste0("fitted", 1:sum(Q != 0))] |> unlist()
  gr.best <- ff$gr(pars.best)
  obj.best <- ff$fn(pars.best)
  fn.cdf <- ecdf(result_good$objective)
  
  result <- list(pars.best = pars.best,
                 obj.best = obj.best,
                 gr.best = gr.best,
                 Phylo = Phylodata,
                 pars.start = pars.start,
                 fn.cdf = fn.cdf)
  
  cat('loglik:\n'); print(obj.best)
  cat('pars:\n'); print(pars.best)
  cat('loglik cdf:\n'); plot(fn.cdf)
  class(result) <- "PAD"
  return(result)
}

# without parallel
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