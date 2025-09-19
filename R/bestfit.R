library(corHMM)
library(Matrix)
library(ape)
library(expm)
library(RTMB)

source('R/Q_template.R')
source('R/getinfo.R')

bestfit  <- function(dat, ..., rate.cat = 1, multistart = 10, jitter.sd = 0.25) {
  tt <- system.time(
    if (multistart > 1) {
      Q <- Q_template(k= ncol(dat$data)-1)
      nrates <- sum(Q != 0)
      p <- rexp(nrates, rate = 1)
      dat <- list(tree = dat$tree, data = dat$data, rate.mat = Q, p = p)
      
      invisible(capture.output(suppressWarnings(x <- with(dat, corHMM(tree, data, rate.mat, p = p, rate.cat = rate.cat, ...)))))
      jitter <- pmin(pmax(rnorm(multistart-1, mean = 0, sd = jitter.sd), 0), 1)
      loglik <- list()
      model <- list()
      loglik[1] <- x$loglik
      model[[1]] <- x
      
      for (i in 1:(multistart-1)) {
        p <- log(jitter[i] + exp(p))
        dat <- list(tree = dat$tree, data = dat$data, rate.mat = Q, p = p)
        invisible(capture.output(suppressWarnings(y <- with(dat, corHMM(tree, data, rate.mat, p = p, rate.cat = rate.cat, ...)))))
        model[[i+1]] <- y
        loglik[i+1] <- y$loglik
      }
      
      best <- model[[which.max(loglik)]]
    })
  attr(best, "time") <- tt
  best
}

out <- simfun(ntrait = 2, ntaxa = 20, seed = 2)
o1 <- bestfit(out)
attr(o1, 'time')
