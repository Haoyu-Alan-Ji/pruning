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