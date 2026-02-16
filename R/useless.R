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