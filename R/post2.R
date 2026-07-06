postfun_tmbstan <- function(pars, Phylodata) {

  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  "diag<-" <- RTMB::ADoverload("diag<-")

  getAll(pars, Phylodata)

  nll <- prune_nll(pars, Phylodata)
  loglik <- -nll

  logdnorm <- function(x, mu, sd) {
    -0.5 * log(2 * pi * sd^2) - 0.5 * ((x - mu) / sd)^2
  }

  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)

  log.prior <- q_par[1] * 0

  if (q_prep$mode == "rates") {
    log.prior <- sum(logdnorm(q_par, prior.mean, prior.sd))
  } else {
    for (nm in names(q_prep$blocks)) {
      b <- q_prep$blocks[[nm]]
      beta <- q_par[b$par_index]
      int <- colnames(b$X) == "(Intercept)"

      if (any(int)) {
        log.prior <- log.prior + sum(logdnorm(beta[int], prior.mean, prior.sd))
      }

      if (any(!int)) {
        log.prior <- log.prior + sum(logdnorm(beta[!int], 0, coef_sd))
      }
    }
  }

  gl.log.prior <- q_par[1] * 0

  if (gainloss != 0 && length(gainloss_pairs) > 0) {
    gl.mean <- (lb_gainloss + ub_gainloss) / 2
    gl.sd <- (ub_gainloss - lb_gainloss) / (2 * range_gainloss)

    gl.values <- q_par[1] * 0 + numeric(length(gainloss_pairs))

    for (i in seq_along(gainloss_pairs)) {
      idx <- gainloss_pairs[[i]]
      gl.values[i] <- q_par[idx[1]] - q_par[idx[2]]
    }

    gl.log.prior <- sum(logdnorm(gl.values, gl.mean, gl.sd))
  }

  logpost <- loglik + log.prior + gainloss * gl.log.prior

  -logpost
}

postMCMC_tmbstan <- function(tree, traitMatrix, formula_list = NULL,
                             state = 2,
                             pars = NULL,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             chains = 4,
                             seed = 427,
                             lb = log(1e-9),
                             ub = log(1e2),
                             range = 3,
                             coef_sd = 2,
                             gainloss = 0,
                             lb_gainloss = log(1e-3),
                             ub_gainloss = log(1e3),
                             range_gainloss = 3,
                             adapt_delta = 0.9,
                             max_treedepth = 12,
                             optimize_first = TRUE) {

  mode <- if (is.null(formula_list)) "rates" else "formula"

  tree <- ape::reorder.phylo(tree, "pruningwise")

  sp <- as.character(traitMatrix[[1]])
  trait_df <- traitMatrix[match(tree$tip.label, sp), -1, drop = FALSE]
  trait_df[] <- lapply(trait_df, as.numeric)

  ntrait <- ncol(trait_df)
  nvec <- if (length(state) == 1) rep(state, ntrait) else state

  w <- rev(cumprod(rev(nvec)))
  w <- c(w[-1], 1)
  trait_values <- rowSums(sweep(as.matrix(trait_df), 2, w, "*")) + 1

  q_prep <- Q_prep(
    mode = mode,
    formula_list = formula_list,
    nstate = state,
    ntrait = ntrait
  )

  if (is.null(pars)) {
    pars <- setNames(rep(log(0.1), q_prep$n_par), q_prep$par_names)
  }

  gainloss_pairs <- gl_pairs(q_prep)

  if (length(gainloss_pairs) == 0) {
    gainloss <- 0
  }

  Phylodata <- list(
    q_prep = q_prep,
    tree = tree,
    trait_values = trait_values,
    gainloss_pairs = gainloss_pairs,

    lb = lb,
    ub = ub,
    range = range,
    coef_sd = coef_sd,

    gainloss = gainloss,
    lb_gainloss = lb_gainloss,
    ub_gainloss = ub_gainloss,
    range_gainloss = range_gainloss
  )

  obj <- RTMB::MakeADFun(
    func = cmb(postfun_tmbstan, Phylodata),
    parameters = list(q_par = pars),
    silent = TRUE
  )

  if (optimize_first) {
      with(obj, nlminb(par, fn, gr))
      cat("fitted pars:\n")
      print(obj$env$last.par.best)
  }
  fit <- tmbstan::tmbstan(
    obj,
    chains = chains,
    iter = iter_warmup + iter_sampling,
    warmup = iter_warmup,
    seed = seed,
    init = if (optimize_first) "last.par.best" else "0",
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )

  list(
    fit = fit,
    obj = obj,
    q_prep = q_prep,
    par_names = q_prep$par_names,
    Phylodata = Phylodata
  )
}
