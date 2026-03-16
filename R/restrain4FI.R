# single formula extraction
if (!inherits(f, "formula")) {
  stop("Each element must be a formula, e.g. ag ~ pc*sc.")
}
if (length(f) != 3L) {
  stop("Each formula must have both a left- and right-hand side.")
}
if (!is.symbol(lhs)) { # only ag, things like log(ag) is not allowed
  stop("Left-hand side must be a single trait name.")
}

# formula list extraction
if (!is.list(formula_list) || length(formula_list) == 0L) {
  stop("formula_list must be a non-empty list of formulas.")
}

if (anyDuplicated(lhs_names)) {
  dup <- unique(lhs_names[duplicated(lhs_names)])
  stop("Duplicate trait formulas found for: ", paste(dup, collapse = ", "))
}

# number of states
if (is.null(names(n))) {
  if (length(n) != length(trait_order)) {
    stop("Unnamed 'n' must have the same length as trait_order.")
  }
  n <- as.integer(n)
  names(n) <- trait_order
  return(n)
}

if (!all(trait_order %in% names(n))) {
  stop("Named 'n' must contain all traits in trait_order.")
}

nstate <- as.integer(nstate[trait_names]) # reorder when c(2, 2, 3) input, for name order is necessary in this situation
names(nstate) <- trait_names

edge_table <- function(traitM, trait_names) {
  ns <- nrow(traitM)
  edge_rows <- vector("list", length = 0L)
  k <- 1L

  Xstate <- as.matrix(traitM[trait_names])

  for (i in seq_len(ns)) {
    for (j in seq_len(ns)) {
      if (i == j) next

      diff_idx <- which(Xstate[i, ] != Xstate[j, ])

      if (length(diff_idx) == 1L) {
        target <- trait_names[diff_idx]
        delta <- Xstate[j, diff_idx] - Xstate[i, diff_idx] #allow 0 to 2 if 3 states

        direction <- if (delta > 0) "gain" else "loss"

        edge_rows[[k]] <- c(
          list(
            edge_id = k, from = i, to = j,
            from_label = traitM$label[i], to_label = traitM$label[j],
            changed_trait = target, direction = direction,
            focal_from = Xstate[i, diff_idx], focal_to = Xstate[j, diff_idx],
            delta = delta
          ),
          as.list(traitM[i, trait_names, drop = FALSE])
        )
        k <- k + 1L
      }
    }
  }
  edge_table <- do.call(rbind, lapply(edge_rows, \(x) as.data.frame(x)))

  int_cols <- c("edge_id", "from", "to", "focal_from", "focal_to", "delta", trait_names)
  edge_table[int_cols] <- lapply(edge_table[int_cols], as.integer)
  edge_table
}

check_formula_variables <- function(formula_list, trait_names) {
  for (tr in names(formula_list)) {
    f <- formula_list[[tr]]
    rhs_vars <- all.vars(f[[3]])

    bad_vars <- setdiff(rhs_vars, trait_names)
    if (length(bad_vars) > 0L) {
      stop(
        "Unknown variable(s) on RHS of ", tr, ": ",
        paste(bad_vars, collapse = ", ")
      )
    }

    if (tr %in% rhs_vars) {
      stop(
        "RHS of ", tr, " must not include the focal trait itself (", tr, ")."
      )
    }
  }
}