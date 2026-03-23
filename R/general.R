#' @examples
#' formula_test <- translate(list(ag ~ pc * sc, pc ~ 1, sc ~ 1))
translate <- function(formula_list, nstate = 2) {
  trait_names <- vapply(formula_list, \(f) as.character(f[[2]]), character(1))
  names(formula_list) <- trait_names

  stateList <- rep.int(as.integer(nstate), length(trait_names))
  names(stateList) <- trait_names

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
  edge_tab <- edge_tab[order(edge_tab$to, edge_tab$from), ]
  edge_tab$edge_id <- seq_len(nrow(edge_tab))

  blocks <- list()
  par_counter <- 1L
  
  edge_groups <- split(edge_tab, list(edge_tab$changed_trait, edge_tab$direction))
  edge_groups <- edge_groups[sort(names(edge_groups))]

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
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- edge_tab$edge_id
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  Q0_sparse <- Matrix::Matrix(0, ns, ns, sparse = TRUE, dimnames = list(traitMatrix$label, traitMatrix$label))

  ## TODO: uncomment this and see what breaks (what's the error?)
  ## Q0 <- matrix(0, ns, ns, dimnames = list(traitMatrix$label, traitMatrix$label))
  ## Q0 <- RTMB::AD(Q0)  ## convert base-R to RTMB/AD type

  list(trait_names = trait_names, stateList = stateList, formulas = formula_list, traitMatrix = traitMatrix,
    edge_table = edge_tab, blocks = blocks, Q_indicator = Q_indicator, Q0 = Q0_sparse,
    n_par = par_counter - 1L, par_names = unlist(lapply(blocks, `[[`, "coef_names"), use.names = FALSE), par_index = lapply(blocks, `[[`, "par_index")
  )
}

#https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
cmb <- function(f, d) function(p) f(p, d)

Q_prep <- function(mode = c("rates", "formula"),
                   formula_list = NULL,
                   nstate = 2,
                   ntrait = 3) {
  mode <- match.arg(mode)

  if (mode == "formula") {
    trans <- translate(formula_list = formula_list, nstate = nstate)

    return(list(mode = "formula", d = nrow(trans$traitMatrix),
                Q0 = trans$Q0, Q_indicator = trans$Q_indicator,
                n_par = trans$n_par, par_names = trans$par_names,
                blocks = trans$blocks, trans = trans
              ))
  }

  ## rates mode: build Q_indicator internally
  Q_indicator <- Q_template(n = nstate, k = ntrait, set_indices = TRUE)
  ns <- nrow(Q_indicator)
  Q0_sparse <- Matrix(0, ns, ns, sparse = TRUE, dimnames = dimnames(Q_indicator))
  list(mode = "rates", d = nrow(Q_indicator),
      Q_indicator = Q_indicator, Q0 = Q0_sparse,
      n_par = max(Q_indicator), par_names = paste0("log_rate_", seq_len(max(Q_indicator)))
  )
}

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

build_Q <- function(q_par, q_prep) {
  Q <- RTMB::AD(q_prep$Q0)
  if (q_prep$mode == "formula") {
    for (nm in names(q_prep$blocks)) {
      b <- q_prep$blocks[[nm]]

      beta_block <- q_par[b$par_index]
      eta_block  <- drop(b$X %*% beta_block)
      rate_block <- exp(eta_block)

      Q[cbind(b$from, b$to)] <- rate_block
    }

  } else if (q_prep$mode == "rates") {
    idx <- which(q_prep$Q_indicator != 0)
    map <- as.integer(q_prep$Q_indicator[idx])
    rates <- exp(q_par)
    Q[idx] <- rates[map]
  } 
  diag(Q) <- -rowSums(Q)
  Q
}

prune_nll <- function(pars, Phylodata) {
  if (!require("RTMB")) stop("install RTMB package")
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)

  ntips <- ape::Ntip(tree)
  if (length(trait_values) != ntips) {
    stop("length(trait_values) must equal the number of tips in tree.")
  }

  d <- q_prep$d
  Q <- build_Q(q_par, q_prep)

  liks <- matrix(NA_real_, nrow = ntips + tree$Nnode, ncol = d)

  ## initialize tip likelihoods
  liks[seq_len(ntips), ] <- 0
  for (i in seq_len(ntips)) {
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
      u <- drop(P %*% childlik) ## FIXME: do we really need this? (remove drop() and see if anything breaks)
      v <- Matrix::drop(v * u)
    }
    comp[i] <- sum(v)
    liks[i, ] <- v / comp[i]
  }

  TIPS <- 1:ape::Ntip(tree)
  root <- ape::Ntip(tree) + 1
  root.p <- rep(1/d, d)
  neg_loglik <- -1*(sum(log(comp[-TIPS])) +log(sum(root.p * liks[root, ])))
  return(neg_loglik)
}
