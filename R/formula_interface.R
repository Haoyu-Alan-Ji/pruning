#' @examples
#' formula_test <- translate(list(ag ~ pc * sc, pc ~ 1, sc ~ 1))
#' formula_test
translate <- function(formula_list, nstate = 2) {
  trait_names <- vapply(formula_list, \(f) as.character(f[[2]]), character(1)); names(formula_list) <- trait_names

  stateList <- rep.int(as.integer(nstate), length(trait_names)); names(stateList) <- trait_names

  traitMatrix <- do.call(expand.grid, lapply(stateList, function(x) 0:(x - 1)))
  traitMatrix$label <- do.call(paste, c(traitMatrix, sep = "|"))

  edge_table <- function(traitM, trait_names) {
    ns <- nrow(traitM)
    Xstate <- as.matrix(traitM[trait_names])
    
    idx_grid <- expand.grid(to = seq_len(ns), from = seq_len(ns))
    idx_grid <- idx_grid[idx_grid$from != idx_grid$to, ]
    
    diff_matrix <- Xstate[idx_grid$to, ] - Xstate[idx_grid$from, ]
    diff_count <- rowSums(diff_matrix != 0)
    
    valid_idx <- which(diff_count == 1L)
    edges <- idx_grid[valid_idx, ]
    final_diffs <- diff_matrix[valid_idx, ]
    
    change_col_idx <- max.col(final_diffs != 0)
    
    res <- data.frame(edge_id = seq_along(valid_idx),
                      from = edges$from, to = edges$to,
                      from_label = traitM$label[edges$from], to_label = traitM$label[edges$to],
                      changed_trait = trait_names[change_col_idx],
                      delta = final_diffs[cbind(seq_along(change_col_idx), change_col_idx)]
                      )
    
    res$direction <- ifelse(res$delta > 0, "gain", "loss")
    res$focal_from <- Xstate[cbind(edges$from, change_col_idx)]
    res$focal_to   <- Xstate[cbind(edges$to, change_col_idx)]
    res <- cbind(res, traitM[edges$from, trait_names])
    
    return(res)
  }
  edge_tab <- edge_table(traitMatrix, trait_names)

  blocks <- list()
  par_counter <- 1L
  
  edge_groups <- split(edge_tab, list(edge_tab$changed_trait, edge_tab$direction))

  for (group_name in names(edge_groups)) {
    dat <- edge_groups[[group_name]]
    if (nrow(dat) == 0) next
    
    tr <- as.character(dat$changed_trait[1])
    dir <- as.character(dat$direction[1])
    
    curr_formula <- formula_list[[tr]]
    rhs_terms <- delete.response(terms(curr_formula))
    
    X <- model.matrix(rhs_terms, data = dat)
    n_col <- ncol(X)
    
    block_name <- paste(tr, dir, sep = "_")
    par_index <- seq(par_counter, length.out = n_col)
    coef_names <- paste0(block_name, ":", colnames(X))

    blocks[[block_name]] <- list(block_name = block_name, trait = tr, direction = dir,
                                  formula = curr_formula, rhs_terms = rhs_terms, edge_id = dat$edge_id,
                                  from = dat$from, to = dat$to,
                                  from_label = dat$from_label, to_label   = dat$to_label,
                                  X = X, n_par = n_col, coef_names = coef_names, par_index = par_index
                                )
    par_counter <- par_counter + n_col
  }

  ns <- nrow(traitMatrix)
  Q_indicator <- matrix(0L, ns, ns)
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- seq_len(nrow(edge_tab))
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  list(trait_names = trait_names, stateList = stateList, formulas = formula_list, traitMatrix = traitMatrix,
    edge_table = edge_tab, blocks = blocks, Q_indicator = Q_indicator,
    n_par = par_counter - 1L, par_names = unlist(lapply(blocks, `[[`, "coef_names"), use.names = FALSE), par_index = lapply(blocks, `[[`, "par_index")
  )
}

#' @examples
#' formula_test <- translate(list(ag ~ pc*sc, pc ~ 1, sc ~ 1))
#' v <- c(
#' # ag_gain: 2 -> 6 (x3) -> 10 (x5) -> 60 (x2)
#' "ag_gain:(Intercept)" = log(2), "ag_gain:pc" = log(3), "ag_gain:sc" = log(5), "ag_gain:pc:sc" = log(2),
#' # ag_loss: 7 -> 14 (x2) -> 21 (x3) -> 84 (x2)
#' "ag_loss:(Intercept)" = log(7), "ag_loss:pc" = log(2), "ag_loss:sc" = log(3), "ag_loss:pc:sc" = log(2),
#' # pc/sc gain/loss constant
#' "pc_gain:(Intercept)" = log(11), "pc_loss:(Intercept)" = log(12),
#' "sc_gain:(Intercept)" = log(13), "sc_loss:(Intercept)" = log(17)
#' )
#' theta_test <- setNames(numeric(formula_test$n_par), formula_test$par_names)
#' theta_test[names(v)] <- v
#' Q <- build_Q(theta_test, formula_test)
#' Q
build_Q <- function(theta, trans_output) {
  if (length(theta) != trans_output$n_par) {
    stop("Length of theta (", length(theta), ") does not match trans_output$n_par (", trans_output$n_par, ")." )
  }

  ns <- nrow(trans_output$traitMatrix)
  labs <- trans_output$traitMatrix$label

  Q0 <- matrix(0, ns, ns, dimnames = list(labs, labs))
  Q <- RTMB::AD(Q0)

  for (nm in names(trans_output$blocks)) {
    b <- trans_output$blocks[[nm]]

    beta_block <- theta[b$par_index]
    eta_block  <- as.vector(b$X %*% beta_block)
    rate_block <- exp(eta_block)

    idx <- cbind(b$from, b$to)
    Q[idx] <- rate_block
  }

  diag(Q) <- -rowSums(Q)
  Q
}

to_bin <- function(x, n = 2) {
  intToBits(x) |> rev() |> as.integer() |> tail(n)
}

#https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
cmb <- function(f, d) function(p) f(p, d)

set.seed(427)
m3 <- ape::rtree(20)
g1 <- reorder(m3, "pruningwise")
pdat <- list(translate_F = formula_test, tree= g1)

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
  d <- nrow(translate_F$traitMatrix)
  repeat {
      s <- phangorn::simSeq(tree, l = 1, Q = formula_test$Q_indicator,
                            type = "USER", levels = seq_len(d), rate = 1)
      if (nrow(unique(as.character(s))) == d) break
    }
  trait_values <- as.numeric(unlist(s))

  Q <- build_Q(theta, translate_F)
  
  liks <- matrix(NA, nrow = ape::Ntip(tree) + tree$Nnode, ncol = d)
  ntips <- length(trait_values)
  ## FIXME: this can be generalized (for polymorphic/unknown states),
  ##  should possibly be moved up a level
  if (nrow(liks) < ntips) {
    stop("Number of columns in 'liks' is less than number of tips (ntips)")
  }
  liks[1:ntips,] <- 0
  for (i in 1:ntips) {
    liks[i, trait_values[i]] <- 1
  }
  comp <- numeric(nrow(liks))
  anc <- unique(tree$edge[, 1])
  for (i in anc) {
    desRows <- which(tree$edge[, 1] == i)
    desNodes <- tree$edge[desRows, 2]
    v <- rep(1, d)
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
  TIPS <- 1:ape::Ntip(tree)
  root <- ape::Ntip(tree) + 1
  root.p <- rep(1/d, d)  
  neg_loglik <- -1*(sum(log(comp[-TIPS])) + log(sum(root.p * liks[root, ])))
  return(neg_loglik)
}

t <- prune_nll(theta_test, pdat)
