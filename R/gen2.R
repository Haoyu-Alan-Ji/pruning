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

Q_prep <- function(mode = c("rates", "formula"), formula_list = NULL, nstate = 2, ntrait = 3) {
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
  Q0_sparse <- Matrix::Matrix(0, ns, ns, sparse = TRUE, dimnames = dimnames(Q_indicator))
  list(mode = "rates", d = nrow(Q_indicator),
      Q_indicator = Q_indicator, Q0 = Q0_sparse,
      n_par = max(Q_indicator), par_names = paste0("log_rate_", seq_len(max(Q_indicator)))
  )
}

QAD <- function(q_prep) {
  d <- q_prep$d

  ## All allowed off-diagonal transitions.
  off_pos <- which(q_prep$Q_indicator != 0L, arr.ind = TRUE)
  off_from <- off_pos[, 1L]
  off_to <- off_pos[, 2L]

  ## Diagonal positions.
  diag_pos <- cbind(seq_len(d), seq_len(d))

  ## Sparse skeleton with all structural entries stored.
  ## Values are placeholders and will be overwritten in build_Q().
  Q_skel <- Matrix::sparseMatrix(
    i = c(off_from, seq_len(d)),
    j = c(off_to, seq_len(d)),
    x = rep(1, length(off_from) + d),
    dims = c(d, d),
    dimnames = dimnames(q_prep$Q_indicator)
  )

  ## Map sparse slot positions.
  ## Matrix stores sparse matrices in compressed-column form.
  slot_i <- Q_skel@i + 1L
  slot_j <- rep(seq_len(ncol(Q_skel)), diff(Q_skel@p))
  slot_key <- paste(slot_i, slot_j, sep = ":")

  off_key <- paste(off_from, off_to, sep = ":")
  diag_key <- paste(seq_len(d), seq_len(d), sep = ":")

  q_prep$Q0_skel <- Q_skel
  q_prep$Q0_AD <- RTMB::AD(Q_skel)

  q_prep$off_pos <- off_pos
  q_prep$off_from <- off_from
  q_prep$off_to <- off_to
  q_prep$off_x_index <- match(off_key, slot_key)
  q_prep$diag_x_index <- match(diag_key, slot_key)

  if (anyNA(q_prep$off_x_index)) {
    stop("Failed to match some off-diagonal entries to sparse @x positions.")
  }

  if (anyNA(q_prep$diag_x_index)) {
    stop("Failed to match some diagonal entries to sparse @x positions.")
  }

  ## For rates mode.
  q_prep$rate_map <- as.integer(q_prep$Q_indicator[off_pos])

  ## For formula mode.
  if (!is.null(q_prep$blocks)) {
    for (nm in names(q_prep$blocks)) {
      b <- q_prep$blocks[[nm]]

      b_key <- paste(b$from, b$to, sep = ":")
      b$x_index <- match(b_key, slot_key)

      if (anyNA(b$x_index)) {
        stop("Failed to match sparse @x positions for block: ", nm)
      }

      q_prep$blocks[[nm]] <- b
    }
  }

  q_prep
}



build_Q <- function(q_par, q_prep) {

  Q <- q_prep$Q0_AD

  ## AD zero with the correct type.
  zero <- q_par[1L] * 0

  ## Start from a zeroed sparse value vector.
  ## In MakeADFun this becomes an advector.
  xval <- Q@x * zero

  if (q_prep$mode == "formula") {
    for (nm in names(q_prep$blocks)) {
      b <- q_prep$blocks[[nm]]

      beta_block <- q_par[b$par_index]
      eta_block <- drop(b$X %*% beta_block)
      rate_block <- exp(eta_block)

      xval[b$x_index] <- rate_block
    }

  } else if (q_prep$mode == "rates") {
    rates <- exp(q_par)
    xval[q_prep$off_x_index] <- rates[q_prep$rate_map]

  } else {
    stop("Unknown `q_prep$mode`: ", q_prep$mode)
  }

  ## Fill diagonal entries as negative row sums of off-diagonal rates.
  row_total <- rep(zero, q_prep$d)

  for (kk in seq_along(q_prep$off_x_index)) {
    rr <- q_prep$off_from[kk]
    row_total[rr] <- row_total[rr] + xval[q_prep$off_x_index[kk]]
  }

  xval[q_prep$diag_x_index] <- -row_total

  ## Critical part:
  ## Ordinary Q@x <- xval fails because dgCMatrix requires numeric @x.
  ## We intentionally bypass S4 validity checking here so that RTMB::expAv()
  ## can see AD-valued sparse entries.
  methods::slot(Q, "x", check = FALSE) <- xval

  Q
}

sparseAD <- function(Q, a) {
  A <- Q

  ## Do not use Q * a.
  ## Sparse Matrix scalar multiplication tries to replace @x
  ## through the ordinary S4 validity path, which rejects advector.
  methods::slot(A, "x", check = FALSE) <- Q@x * a

  A
}

prune_nll <- function(pars, Phylodata) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  "diag<-" <- RTMB::ADoverload("diag<-")

  RTMB::getAll(pars, Phylodata)

  ntips <- ape::Ntip(tree)

  if (length(trait_values) != ntips) {
    stop("length(trait_values) must equal the number of tips in tree.")
  }

  d <- q_prep$d
  Q <- build_Q(q_par, q_prep)

  liks <- matrix(
    NA_real_,
    nrow = ntips + tree$Nnode,
    ncol = d
  )

  ## Initialize tip likelihoods.
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

      ## Apply exp(Q * t) to the child likelihood without forming P.
      ## Requires Q to be a sparse Matrix-like object with @x.
      Qt <- sparseAD(Q, t)
      u <- RTMB::expAv(Qt, as.vector(childlik))

      ## Keep your current shape-control logic.
      v <- drop(as.matrix(v * u))
    }

    comp[i] <- sum(v)
    liks[i, ] <- v / comp[i]
  }

  TIPS <- seq_len(ape::Ntip(tree))
  root <- ape::Ntip(tree) + 1L
  root.p <- rep(1 / d, d)

  neg_loglik <- -1 * (
    sum(log(comp[-TIPS])) +
      log(sum(root.p * liks[root, ]))
  )

  neg_loglik
}


if (FALSE) {
  library(ape)
  library(Matrix)
  library(RTMB)

  set.seed(1)

  check_fit <- function(q_prep, tree, label = "test") {
    q_prep <- QAD(q_prep)

    trait_values <- sample(
      seq_len(q_prep$d),
      size = ape::Ntip(tree),
      replace = TRUE
    )

    Phylodata <- list(
      q_prep = q_prep,
      tree = tree,
      trait_values = trait_values
    )

    pars <- list(
      q_par = setNames(
        rep(log(0.1), q_prep$n_par),
        q_prep$par_names
      )
    )

    cat("\n---", label, ": direct build_Q ---\n")
    Q0 <- build_Q(pars$q_par, q_prep)
    print(Q0)
    stopifnot(nrow(Q0) == q_prep$d)
    stopifnot(ncol(Q0) == q_prep$d)

    cat("\n---", label, ": direct prune_nll ---\n")
    nll0 <- prune_nll(pars, Phylodata)
    print(nll0)
    stopifnot(is.finite(as.numeric(nll0)))

    cat("\n---", label, ": MakeADFun ---\n")
    ff <- RTMB::MakeADFun(
      func = cmb(prune_nll, Phylodata),
      parameters = pars,
      silent = TRUE
    )

    cat("\n---", label, ": fn ---\n")
    fn0 <- ff$fn(ff$par)
    print(fn0)

    if (!is.finite(fn0)) {
      print(warnings())
      stop("Non-finite objective from ff$fn(ff$par).")
    }

    cat("\n---", label, ": gr ---\n")
    gr0 <- ff$gr(ff$par)
    print(gr0)
    stopifnot(length(gr0) == q_prep$n_par)
    stopifnot(all(is.finite(gr0)))

    cat("\n---", label, ": short nlminb smoke test ---\n")
    opt <- nlminb(
      start = ff$par,
      objective = ff$fn,
      gradient = ff$gr,
      control = list(iter.max = 5, eval.max = 10)
    )

    print(opt$convergence)
    print(opt$objective)

    invisible(list(
      q_prep = q_prep,
      Phylodata = Phylodata,
      ff = ff,
      opt = opt
    ))
  }

  tree <- ape::rtree(10)
  tree <- ape::reorder.phylo(tree, "pruningwise")

  formula_list <- list(
    A ~ B + C,
    B ~ A + C,
    C ~ A + B
  )

  q_prep_formula <- Q_prep(
    mode = "formula",
    formula_list = formula_list,
    nstate = 2,
    ntrait = 3
  )

  fit_formula_test <- check_fit(
    q_prep = q_prep_formula,
    tree = tree,
    label = "formula mode"
  )

  q_prep_rates <- Q_prep(
    mode = "rates",
    nstate = 2,
    ntrait = 3
  )

  fit_rates_test <- check_fit(
    q_prep = q_prep_rates,
    tree = tree,
    label = "rates mode"
  )

  message("\nAll moved-AD-template smoke tests passed.")
}