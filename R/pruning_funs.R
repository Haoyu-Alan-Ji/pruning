#' @param trans_rates vector of transition rates
#' @param Q_template a square matrix with zero values on the diagonal, and corresponding to any forbidden transitions, and integer indices in positions corresponding to distinct transition rates
setup_Q <- function(trans_rates, Q_template) {
  Q <- Q_template
  Q[Q!=0] <- trans_rates[Q[Q!=0]] 
  diag(Q) <- -1*rowSums(Q)
  return(Q)
}

#' @param pars named list of parameters (
#' TMBdata should be a named list containing `tree` (a `phylo` object as defined in the `ape` package);
#' trait_values (a vector of integer trait values, corresponding to states at tips
prune_nll <- function(pars) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, TMBdata)
  Q <- setup_Q(exp(log_trans_rates), Q_template)
  liks <- matrix(NA, nrow = Ntip(tree) + tree$Nnode, ncol = d)
  ntips <- length(trait_values)
  ## FIXME: this can be generalized (for polymorphic/unknown states),
  ##  should possibly be moved up a level
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

hmm2prune <- function(corHMM_fit){
  index.mat <- corHMM_fit$index.mat
  solution <- corHMM_fit$solution
  
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

#' could do this faster (probably) with an iterated Kronecker product, but it's totally unnecessary
#' for any practical use case we have ...
#' @param n (vector of) number of states (if n is a scalar, all traits have the same number of states)
#' @param k number of traits (if n is a vector, k will be `length(n)`
setup_Q_template <- function(n=3, k= 1) {
  if (length(n) == 1) {
    n <- rep(n, k)
  }
  all_states <- do.call(expand.grid, lapply(n, \(x) 0:(x-1)))
  ns <- prod(n)
  m <- matrix(0, ns, ns)
  for (i in 1:ns) {
    ## exactly one state changes ...
    for (j in 1:ns) {
      m[i,j] <- as.numeric(sum(all_states[i,] != all_states[j,])== 1)
    }
  }
  return(m)
}

  
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

if (FALSE) {
  ## testing
  system.time(Q_big <- setup_Q_template(k=5)) ## 10 seconds
  png("bigmat.png")
  imat(Q_big)
  dev.off()
}

