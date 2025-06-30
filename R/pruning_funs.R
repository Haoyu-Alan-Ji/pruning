#' @param n (vector of) number of states (if n is a scalar, all traits have the same number of states)
#' @param k number of traits (only used if n is a scalar)
#' @examples
#' setup_Q_template(n = 3, k = 2)
#' setup_Q_template(n = 2, k = 3)
setup_Q_template <- function(n=3, k= 1, set_indices = TRUE) {
  if (length(n) == 1) {
    n <- rep(n, k)
  }
  all_states <- do.call(expand.grid, lapply(n, \(x) 0:(x-1)))
  dimnms <- apply(all_states, 1, \(x) sprintf("(%s)", paste(x, collapse = ",")))
  ns <- prod(n)
  m <- matrix(0, ns, ns)
  for (i in 1:ns) {
    ## exactly one state changes ...
    for (j in 1:ns) {
      m[i,j] <- as.numeric(sum(all_states[i,] != all_states[j,])== 1)
    }
  }
  dimnames(m) <- list(dimnms, dimnms)
  if (set_indices) {
    m[m!=0] <- seq_len(sum(m==1))
  }
  return(m)
}

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

#' @param trans_rates vector of transition rates
#' @param Q_template a square matrix with zero values on the diagonal, and corresponding to any forbidden transitions, and integer indices in positions corresponding to distinct transition rates
setup_Q <- function(trans_rates, Q_template) {
  Q <- Q_template
  Q[Q!=0] <- trans_rates[Q[Q!=0]] 
  diag(Q) <- -1*rowSums(Q)
  return(Q)
}

#' @param pars named list of parameters (
#' TMBdata *should exist in the environment* and should
#' be a named list containing
#' - `tree` (a `phylo` object as defined in the `ape` package)
#' - trait_values (a vector of integer trait values, corresponding to states at tips)
#' - Q_template (a matrix with non-zero integer indices at all allowed locations)
prune_nll <- function(pars) {
  if (!require("RTMB")) stop("install RTMB package")
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, TMBdata)
  d <- ncol(Q_template)
  Q <- setup_Q(exp(log_trans_rates), Q_template)
  liks <- matrix(NA, nrow = Ntip(tree) + tree$Nnode, ncol = d)
  ntips <- length(trait_values)
  ## FIXME: this can be generalized (for polymorphic/unknown states),
  ##  should possibly be moved up a level
  if (nrow(liks) < ntips) {
    stop("Number of columns in 'liks' is less than number of tips (ntips)")
  }
  liks[1:ntips,] <- 0
  for (i in 1:nrow(liks)) {
    liks[i, trait_values[i]] <- 1
  }
  comp <- numeric(nrow(liks))
  anc <- unique(tree$edge[, 1])
  for (i in anc) {
    desRows <- which(tree$edge[, 1] == i)
    desNodes <- tree$edge[desRows, 2]
    v <- 1
    for (j in seq_along(desRows)) {
      t <- tree$edge.length[desRows[j]]
      P <- Matrix::expm(Q * t)
      v <- drop(v * (P %*% liks[desNodes[j], ]))
    }
    comp[i] <- sum(v)
    liks[i, ] <- v / comp[i]
  }
  TIPS <- 1:Ntip(tree)
  root <- Ntip(tree) + 1
  root.p <- rep(1/d, d)  
  neg_loglik <- -1*(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root, ])))
  return(neg_loglik)
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
