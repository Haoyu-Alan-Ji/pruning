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

setup_Q_template <- function(ntraits=3, nstates = 2) {
  ## machinery from corHMM/R/rate.mat.maker.R
  ## in corHMM, nl = nstates, nk = ntraits
  ## for now, n binary traits [expand to allow different numbers of states per trait??]
  sn0 <- c(t(outer(0:(n-1), 0:(n-1), paste, sep = ",")))
  states <- sprintf("(%s)", sn0)
  mdim <- ntraits^nstates
  matlist <- replicate(ntraits, matrix(NA, mdim, mdim))
  ## binary digits
  vecs 
  Q_template <- matrix(0, length(states), length(states),
                       dimnames = list(states, states))
  ## how do we programmatically determine allowed states?

  ## is there machinery in the 
  allowed <- matrix(c(2,1,
                    3,1,
                    1,2,
                    4,2,
                    1,3,
                    4,3,
                    3,4,
                    2,4),
                  ncol = 2,
                  byrow = TRUE)
  Q_template[allowed] <- 1:8
}

nl <- 2
k <- 3
mat1 <- matrix(, nl^k, nl^k)
mat2 <- matrix(, nl^k, nl^k)
mat3 <- matrix(, nl^k, nl^k)
      vec.tmp1 <- c(0, 1, 0, 0, 1, 1, 0, 1)
      vec.tmp2 <- c(0, 0, 1, 0, 1, 0, 1, 1)
      vec.tmp3 <- c(0, 0, 0, 1, 0, 1, 1, 1)
      for (i in 1:(nl^k)) {
        mat1[i, ] <- abs(vec.tmp1 - vec.tmp1[i])
        mat2[i, ] <- abs(vec.tmp2 - vec.tmp2[i])
        mat3[i, ] <- abs(vec.tmp3 - vec.tmp3[i])
      }
      matFINAL <- mat1 + mat2 + mat3
