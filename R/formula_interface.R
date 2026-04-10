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

#' take a list of two-sided formulas. Check whether they contain sym() on the LHS;
#' if so, set `is_sym` to TRUE and remove sym()
#' @return a list containing `formula_list` (processed formula list) and `is_sym` (whether each term
#' contains `sym()`
#' @examples
#' get_sym(list(A~sym(a+b*c), B~sym(1), C~a+b))
get_sym <- function(formula_list) {
  is_sym <- rep(FALSE, length(formula_list))
  for (i in seq_along(formula_list)) {
    f <- formula_list[[i]]
    if (identical(f[[3]][[1]], quote(sym))) {
      is_sym[i] <- TRUE
      f[[3]] <- f[[3]][[2]] ## drop
      formula_list[[i]] <- f
    }
  }
  return(list(formula_list = formula_list, is_sym = is_sym))
}

#' @examples
#' formula_test <- translate(list(ag ~ pc * sc, pc ~ 1, sc ~ 1))
#' formula_test2 <- translate(list(ag ~ sym(pc * sc), pc ~ 1, sc ~ 1))
translate <- function(formula_list, nstate = 2) {

  ## process formula list to get symmetry info
  ff <- get_sym(formula_list)
  is_sym <- ff$is_sym
  formula_list <- ff$formula_list
  
  trait_names <- vapply(formula_list, \(f) as.character(f[[2]]), character(1))
  names(formula_list) <- names(is_sym) <- trait_names

  
  ## if symmetric, has one direction been processed?
  sym_done <- rep(FALSE, length(trait_names))
  names(sym_done) <- trait_names
  
  stateList <- rep.int(as.integer(nstate), length(trait_names))
  names(stateList) <- trait_names

  traitMatrix <- do.call(expand.grid, lapply(stateList, function(x) 0:(x - 1)))
  traitMatrix$label <- do.call(paste, c(traitMatrix, sep = "|"))

  edge_tab <- edge_table(traitMatrix, trait_names)
  edge_tab <- edge_tab[order(edge_tab$to, edge_tab$from), ]
  edge_tab$edge_id <- seq_len(nrow(edge_tab))

  blocks <- list()
  par_counter <- 1L
  
  edge_groups <- split(edge_tab, list(edge_tab$changed_trait, edge_tab$direction))
  edge_groups <- edge_groups[sort(names(edge_groups))]

  last_block <- NULL
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

    ## if sym *and* second block for this trait, set par index same as previous;
    ## don't advance par_counter
    if (!sym_done[[tr]]) {
      par_index <- seq(par_counter, length.out = n_col)
      par_counter <- par_counter + n_col
      if (is_sym[[tr]]) sym_done[[tr]] <- TRUE
      last_block <- block_name
    } else {
      if (is_sym[[tr]] && sym_done[[tr]]) {
        browser()
        ## untested!
        par_index <- blocks[[last_block]]$par_index
      }
    }
    coef_names <- paste0(block_name, ":", colnames(X))

    blocks[[block_name]] <- list(block_name = block_name, trait = tr, direction = dir,
                                  formula = curr_formula, rhs_terms = rhs_terms, edge_id = dat$edge_id,
                                  from = dat$from, to = dat$to,
                                  from_label = dat$from_label, to_label   = dat$to_label,
                                  X = X, n_par = n_col, coef_names = coef_names, par_index = par_index
                                )

  }

  ns <- nrow(traitMatrix)
  Q_indicator <- matrix(0L, ns, ns)
  Q_indicator[as.matrix(edge_tab[, c("from", "to")])] <- edge_tab$edge_id
  rownames(Q_indicator) <- colnames(Q_indicator) <- traitMatrix$label

  Q0_sparse <- Matrix::Matrix(0, ns, ns, sparse = TRUE, dimnames = list(traitMatrix$label, traitMatrix$label))
  #Q0 <- matrix(0, ns, ns, dimnames = list(traitMatrix$label, traitMatrix$label))
  #Q0 <- RTMB::AD(Q0)  ## convert base-R to RTMB/AD type

  list(trait_names = trait_names,
       stateList = stateList,
       formulas = formula_list,
       traitMatrix = traitMatrix,
       edge_table = edge_tab,
       blocks = blocks,
       Q_indicator = Q_indicator,
       Q0 = Q0_sparse,
       n_par = par_counter - 1L,
       par_names = unlist(lapply(blocks, `[[`, "coef_names"), use.names = FALSE),
       par_index = lapply(blocks, `[[`, "par_index")
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
## FIXME? should we construct Q once, then insert new values into the matrix
## during the optimization process? (Would still have to do -rowSums(Q))

## FINISHED, now sparse Q construction is operated in translate()
## and AD class giving will do in build_Q() (do this in translate() will lead terrible bug).
## bug description (when give AD in translate()): 
## Q is already AD type, so it's advector, however, Matrix::Matrix() will only operate on Matrix object, thus, error
## In RTMB file, the recommended approach is to "first create a standard object, and then use `AD()` to imbue it with AD semantics",
## rather than feeding a dense matrix (which is *already* an AD object) into `Matrix::Matrix()` to convert it into a sparse matrix.
build_Q <- function(theta, trans_output) {
  if (length(theta) != trans_output$n_par) {
    stop("Length of theta (", length(theta), ") does not match trans_output$n_par (", trans_output$n_par, ")." )
  }
  Q <- RTMB::AD(trans_output$Q0)
  for (nm in names(trans_output$blocks)) {
    b <- trans_output$blocks[[nm]]

    beta_block <- theta[b$par_index]
    eta_block  <- drop(b$X %*% beta_block)
    rate_block <- exp(eta_block)

    idx <- cbind(b$from, b$to)
    Q[idx] <- rate_block
  }
  diag(Q) <- -rowSums(Q)
  Q
}

# to_bin <- function(x, n = 2) {
#   intToBits(x) |> rev() |> as.integer() |> tail(n)
# }

#https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
cmb <- function(f, d) function(p) f(p, d)

#' @param pars named list of parameters (
#' Phylodata *should exist in the environment* and should
#' be a named list containing
#' - `tree` (a `phylo` object as defined in the `ape` package)
#' - trait_values (a vector of integer trait values, corresponding to states at tips)
#' - translated formula (the output contains almost everything we want)
#' @examples
#' set.seed(427)
#' m3 <- ape::rtree(20)
#' g1 <- reorder(m3, "pruningwise")
#' d <- nrow(formula_test$traitMatrix)
#' repeat {
#'     s <- phangorn::simSeq(g1, l = 1, Q = formula_test$Q_indicator,
#'                           type = "USER", levels = seq_len(d), rate = 1)
#'     if (nrow(unique(as.character(s))) == d) break
#'   }
#' s <- as.numeric(unlist(s))
#' pdat <- list(translate_F = formula_test, tree= g1, trait_values = s)
#' testval <- prune_nll(list(theta = theta_test), pdat)
#' ff <- RTMB::MakeADFun( func = cmb(prune_nll, pdat),  parameters = list(theta = theta_test), silent = TRUE)
#' ff$fn(); ff$gr()
#' opt1 <- with(ff, nlminb(par, fn, gr))
## FIXME: rearrange so that simSeq and build_Q (without computing values to insert)
## should be outside of the function that will be passed to MakeADFun()
## i.e. they should be part of Phylodata instead  

## FINISHED, see @param
prune_nll <- function(pars, Phylodata) {
  if (!require("RTMB")) stop("install RTMB package")
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(pars, Phylodata)

  d <- nrow(translate_F$traitMatrix)
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
    v <- 1 # rep(1, d) would also work, but more stable...?
    for (j in seq_along(desRows)) {
      t <- tree$edge.length[desRows[j]]
      childlik <- matrix(liks[desNodes[j], ], ncol = 1)
      P <- Matrix::expm(Q * t)     
      ## FIXME: try u <- RTMB::expAv(Q * t, as.vector(childlik))??
      u <- drop(P %*% childlik)
      ## need as.matrix() to make drop() work properly 
      v <- Matrix::drop(v * u)
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

make_prune_problem <- function(tree, trait_values,
                               mode = c("rates", "formula"),
                               formula_list = NULL,
                               nstate = 2L,
                               Q_template = NULL,
                               start = NULL) {
  mode <- match.arg(mode)

  q_spec <- prepare_q_spec(
    mode = mode,
    formula_list = formula_list,
    nstate = nstate,
    Q_template = Q_template
  )

  par0 <- setNames(numeric(q_spec$n_par), q_spec$par_names)

  if (!is.null(start)) {
    if (is.null(names(start))) {
      if (length(start) != length(par0)) {
        stop("Unnamed start vector has wrong length.")
      }
      par0[] <- start
    } else {
      par0[names(start)] <- start
    }
  }

  pdat <- list(
    tree = tree,
    trait_values = trait_values,
    q_spec = q_spec
  )

  list(
    pdat = pdat,
    par0 = par0,
    q_spec = q_spec
  )
}

make_prune_objective <- function(tree, trait_values,
                                 mode = c("rates", "formula"),
                                 formula_list = NULL,
                                 nstate = 2L,
                                 Q_template = NULL,
                                 start = NULL,
                                 silent = TRUE) {
  prob <- make_prune_problem(
    tree = tree,
    trait_values = trait_values,
    mode = mode,
    formula_list = formula_list,
    nstate = nstate,
    Q_template = Q_template,
    start = start
  )

  obj <- RTMB::MakeADFun(
    func = cmb(prune_nll_general, prob$pdat),
    parameters = list(q_par = prob$par0),
    silent = silent
  )

  list(
    obj = obj,
    pdat = prob$pdat,
    q_spec = prob$q_spec,
    par0 = prob$par0
  )
}

if (FALSE) {
    ## testing
    Qtemp <- Q_template(n = 2, k = 3)

    fit1 <- make_prune_objective(
        tree = g1,
        trait_values = s,
        mode = "rates",
        Q_template = Qtemp,
        silent = TRUE
    )

    fit1$obj$fn()
    fit1$obj$gr()
    opt1 <- with(fit1$obj, nlminb(par, fn, gr))
}
