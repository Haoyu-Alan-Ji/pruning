#' @examples
#' streams <- make_streams(9, seed6 = c(1,2,3,4,5,6))
#' sapply(streams, \(s) rstream.sample(s, 3))
library(rstream)

streams <- function(n, seed6 = c(1,2,3,4,5,6),
                         antithetic = FALSE, incprecision = FALSE) {
  stopifnot(length(seed6) == 6)
  base <- new("rstream.mrg32k3a", seed = seed6,
              force.seed = TRUE, antithetic = antithetic,
              incprecision = incprecision)

  out <- vector("list", n)

  out[[1]] <- rstream.clone(base)
  rstream.nextsubstream(out[[1]])

  if (n > 1) {
    for (i in 2:n) {
      out[[i]] <- rstream.clone(out[[i-1]])
      rstream.nextsubstream(out[[i]])
    }
  }
  out
}

