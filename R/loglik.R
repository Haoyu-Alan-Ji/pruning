#' @param trans_rates vector of transition rates
#' @param Q_template a square matrix with zero values on the diagonal, and corresponding to any forbidden transitions, and integer indices in positions corresponding to distinct transition rates
setup_Q <- function(trans_rates, Q_template) {
  Q <- Q_template
  Q[Q!=0] <- trans_rates[Q[Q!=0]] 
  diag(Q) <- -1*rowSums(Q)
  return(Q)
}

#https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
cmb <- function(f, d) function(p) f(p, d)

#' @param pars named list of parameters (
#' TMBdata *should exist in the environment* and should
#' be a named list containing
#' - `tree` (a `phylo` object as defined in the `ape` package)
#' - trait_values (a vector of integer trait values, corresponding to states at tips)
#' - Q_template (a matrix with non-zero integer indices at all allowed locations)
prune_nll <- function(pars, Phylodata) {
  if (!require("RTMB")) stop("install RTMB package")
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)
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