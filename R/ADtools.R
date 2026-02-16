#' @param n (vector of) number of states (if n is a scalar, all traits have the same number of states)
#' @param k number of traits (only used if n is a scalar)
#' @examples
#' setup_Q_template(n = 3, k = 2)
#' setup_Q_template(n = 2, k = 3)
Q_template <- function(n=2, k= 3, set_indices = TRUE) {
  if (length(n) == 1) {
    n <- rep(n, k)
  }
  all_states <- do.call(expand.grid, lapply(n, \(x) 0:(x-1)))
  ## dimnms to match corHMM standard
  dimnms <- apply(all_states, 1, \(x) paste(x, collapse = "|"))
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

to_bin <- function(x, n = 2) {
  intToBits(x) |> rev() |> as.integer() |> tail(n)
}


# The expm method of Matrix cannot be optimized by AD, but RTMB can accept this, provided that Q for signature 'advector'. https://rdrr.io/cran/RTMB/man/ADmatrix.html

#' @param trans_rates vector of transition rates
#' @param Q_template a square matrix with zero values on the diagonal, and corresponding to any forbidden transitions, and integer indices in positions corresponding to distinct transition rates
setup_Q <- function(log_trans_rates, Q_template) {
  idx <- which(Q_template != 0)      
  map <- as.integer(Q_template[idx])   
  rates <- exp(log_trans_rates)      
  
  Q <- RTMB::AD(matrix(0, nrow(Q_template), ncol(Q_template)))
  Q[idx] <- rates[map]
  diag(Q) <- -rowSums(Q)
  Q
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
  Q <- setup_Q(log_trans_rates, Q_template)
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
      childlik <- matrix(liks[desNodes[j], ], ncol = 1)
      P <- Matrix::expm(Q * t)     
      u <- drop(P %*% childlik)     
      v <- drop(v * u)
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