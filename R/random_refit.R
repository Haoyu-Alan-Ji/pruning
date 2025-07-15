random_refit <- function(ff) {
  n <- length(ff$par)
  start <- log(abs(rnorm(n)))
  res <- with(ff, nlminb(start, fn, gr))
  names(start) <- paste0("start", seq_along(start))
  fitted <-  res$par
  names(fitted) <- paste0("fitted", seq_along(fitted))
  with(res, data.frame(objective, convergence, message,
                       rbind(start), rbind(fitted)))
}
