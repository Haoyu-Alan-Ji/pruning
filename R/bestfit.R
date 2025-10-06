library(corHMM)
library(Matrix)
library(ape)
library(expm)
library(RTMB)

source(here::here('R', 'Q_template.R'))
source(here::here('R', 'getinfo.R'))

bestfit  <- function(dat, ..., rate.cat = 1, multistart = 10, jitter.sd = 0.25) {
  require(c(parallel, rstream))
  tt <- system.time(
    Q <- Q_template(k= ncol(dat$data)-1)
    nrates <- sum(Q != 0)
    p <- rexp(nrates, rate = 1)
    dat <- list(tree = dat$tree, data = dat$data, rate.mat = Q, p = p)

    ## create independent RNG streams
    streams <- lapply(1:(multistart - 1), function(i) rstream.mrg(seed = i + sample.int(1e6, 1)))
    
    ## worker function
    fitfun <- function(i, p_init, dat, jitter.sd, rate.cat, ...) {
      setRNG(streams[[i]])  # unique RNG per worker
      jitter <- pmin(pmax(rnorm(1, mean = 0, sd = jitter.sd), 0), 1)
      p_new <- log(jitter + exp(p_init))
      dat <- list(tree = dat$tree, data = dat$data, rate.mat = dat$rate.mat, p = p_new)
      invisible(capture.output(suppressWarnings(
        y <- with(dat, corHMM(tree, data, rate.mat, p = p, rate.cat = rate.cat, ...))
      )))
      return(list(loglik = y$loglik, model = y))
    }
    
    ## run in parallel
    cl <- makeCluster(min(multistart - 1, detectCores() - 1))
    clusterExport(cl, c("fitfun", "dat0", "p_init", "rate.cat", "jitter.sd", "streams", "..."), 
                  envir = environment())
    parallel_out <- parLapply(cl, 1:(multistart - 1), function(i)
      fitfun(i, p_init, dat0, jitter.sd, rate.cat, ...)
    )
    stopCluster(cl)
    
    ## collect results
    for (i in seq_along(parallel_out)) {
      loglik[i + 1] <- parallel_out[[i]]$loglik
      model[[i + 1]]  <- parallel_out[[i]]$model
    }
    
    best <- model[[which.max(loglik)]]
  )
  attr(best, "time") <- tt
  best
}
