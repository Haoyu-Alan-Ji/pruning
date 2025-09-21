library(corHMM)
library(rbenchmark)

parfun <- function(x) log(na.omit(c(x$solution)))

fitfun  <- function(dat, ..., rate.cat = 1) {
  tt <- system.time(
    invisible(capture.output(suppressWarnings(x <- with(dat, corHMM(tree, data, rate.cat = rate.cat, ...)))))
  )
  attr(x, "time") <- tt
  x
}

bestfit  <- function(dat, ..., rate.cat = 1, multistart = 10, jitter.sd = 0.25) {
  tt <- system.time(
    if (multistart > 1) {
      Q <- Q_template(k= ncol(dat$data)-1)
      nrates <- sum(Q != 0)
      p <- rexp(nrates, rate = 1)
      dat <- list(tree = dat$tree, data = dat$data, rate.mat = Q, p = p)
      
      invisible(capture.output(suppressWarnings(
        x <- with(dat, corHMM(tree, data, rate.mat, p = p, rate.cat = rate.cat, ...)))))
      jitter <- pmin(pmax(rnorm(multistart-1, mean = 0, sd = jitter.sd), 0), 1)
      loglik <- list()
      model <- list()
      loglik[1] <- x$loglik
      model[[1]] <- x
      
      for (i in 1:(multistart-1)) {
        p <- log(jitter[i] + exp(p))
        dat <- list(tree = dat$tree, data = dat$data, rate.mat = Q, p = p)
        invisible(capture.output(suppressWarnings(
          y <- with(dat, corHMM(tree, data, rate.mat, p = p, rate.cat = rate.cat, ...)))))
        model[[i+1]] <- y
        loglik[i+1] <- y$loglik
      }
      
      best <- model[[which.max(loglik)]]
    })
  attr(best, "time") <- tt
  best
}

sumfun <- function(ntrait = 2, ntaxa = 200, model = "ARD", multistart = FALSE, 
                   try.times = 10, jitter.sd = 0.25, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  seed <- seed %||% NA
  ss <- simfun(ntrait = ntrait, ntaxa = ntaxa, seed = seed)
  if (multistart == TRUE) {
    fit_orig <- bestfit(ss, model = model, 
                        multistart = try.times, jitter.sd = jitter.sd, ...)
    fit_RTMB <- bestfit(ss, use_RTMB = TRUE, model = model, 
                        multistart = try.times, jitter.sd = jitter.sd, ...)
  } else {
    fit_orig <- fitfun(ss, model = model, ...)
    fit_RTMB <- fitfun(ss, use_RTMB = TRUE, model = model, ...)
    }
  p_orig <- parfun(fit_orig)
  p_RTMB <- parfun(fit_RTMB)
  p_diff <- p_orig-p_RTMB
  ## truncate par values at -10 (these diffs are mostly irrelevant)
  p_orig_trunc <- pmax(p_orig, -10)
  p_RTMB_trunc <- pmax(p_RTMB, -10)
  p_diff_trunc <- p_orig_trunc - p_RTMB_trunc
  ## get RMSE vs true rates on truncated scale ...
  p_RTMB_rmse <- sqrt(mean((p_RTMB_trunc - ss$true_rates)^2))
  p_orig_rmse <- sqrt(mean((p_orig_trunc - ss$true_rates)^2))
  data.frame(seed, ntrait, ntaxa, model,
             ## RTMB_opt.time = fit_RTMB$opt.time[["elapsed"]],
             ## orig_opt.time = fit_orig$opt.time[["elapsed"]],
             RTMB_tot.time = attr(fit_RTMB, "time")[["elapsed"]],
             orig_tot.time = attr(fit_orig, "time")[["elapsed"]],
             RTMB_loglik = fit_RTMB$loglik,
             orig_loglik = fit_orig$loglik,
             par.rmsdiff = sqrt(mean(p_diff^2)),
             par.maxdiff = max(abs(p_diff)),
             par.rmsdiff.trunc = sqrt(mean(p_diff_trunc^2)),
             par.maxdiff.trunc = max(abs(p_diff_trunc)),
             RTMB_rmse = p_RTMB_rmse,
             orig_rmse = p_orig_rmse
  )
}