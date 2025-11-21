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

## setname <- function(x) {cbind(nm = rownames(x), x)}
##' @param nstate number of states per trait (currently limited to 2)
##' @param ntrait number of traits
##' @param ntaxa number of taxa/phylogenty tips
##' @param seed random-number seed
##' @param meanrate mean transition rate
##' @param seql trait
##' @param collapse_allow for corHMM collapse function
##' @examples
##' simfun(ntaxa = 8)
simfun <- function(nstate = 2, ntrait = 2, ntaxa = 20, seed = NULL,
                   meanrate = 1, seql = 1, collapse = FALSE) {
  if (nstate!=2) stop("oops, simSeq to trait matrix not implemented for nstate!=2, 
                      to_bin's binary algorithm determines that the state can only be 2")
  if (!is.null(seed)) set.seed(seed)
  require("ape")
  require("phangorn")
  phy <- ape::rtree(ntaxa)
  phy <- reorder(phy, "pruningwise")
  
  Q <- Q_template(n=nstate, k= ntrait)
  nrates <- sum(Q != 0)
  Q[Q!=0] <- rexp(nrates, rate = 1/meanrate)
  
  s <- phangorn::simSeq(phy, l = seql, Q = Q,
                        type = "USER", levels = seq(nstate^ntrait),
                        rate = 1)
  if (collapse == FALSE) {
    repeat {
      s <- phangorn::simSeq(phy, l = seql, Q = Q,
                            type = "USER", levels = seq(nstate^ntrait),
                            rate = 1)
      if (nrow(unique(as.character(s))) == prod(nstate^ntrait)) break
    }
  }
  traitMatrix <- sapply(s, function(x) to_bin(x = x-1, n = ntrait)) |>
    t() |>
    as.data.frame() #|>
  #setname()
  
  list(tree = phy, data = traitMatrix)
}

realfun <- function(tree) {
  trait_number = ncol(tree$edge) + 1
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

o1 <- realfun(fish_tree)
