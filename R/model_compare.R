## TO CHECK: check against parameters? against full matrix?
#' @param corHMM_fit a fitted corHMM object

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