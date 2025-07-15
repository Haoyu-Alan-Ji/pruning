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