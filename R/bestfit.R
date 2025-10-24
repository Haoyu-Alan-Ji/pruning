library(corHMM)
library(Matrix)
library(ape)
library(expm)
library(RTMB)
library(future.apply)
library(rstream)
library(parallel)

source(here::here('R', 'Q_template.R'))
source(here::here('R', 'getinfo.R'))

## fun
bestfit  <- function(dat, ..., rate.cat = 1, multistart = 10) {

  tt <- system.time({
    dat <- list(tree = dat$tree, data = dat$data)

    invisible(capture.output(suppressWarnings(
      x <- with(dat, corHMM(tree, data, rate.cat = rate.cat, ...))
    )))
    
    loglik <- list(); loglik[1] <- x$loglik
    model <- list(); model[[1]] <- x

    ## worker function
    fitj <- function(i, dat, rate.cat, ...) {
      s <- new("rstream.mrg32k3a",
           seed = c(i, i+1, i+2, i+3, i+4, i+5),
           force.seed = TRUE)
      rstream.RNG(s)

      cat("▶ Worker", i, "started on PID", Sys.getpid(), "\n")
      flush.console()
      
      dat0 <- list(tree = dat$tree, data = dat$data)
      invisible(capture.output(suppressWarnings(
        y <- with(dat0, corHMM(tree, data, rate.cat = rate.cat, ...))
      )))
      return(list(loglik = y$loglik, model = y))
    }
    
    parallel_out <- future_lapply(
        1:(multistart - 1), \(i) fitj(i, dat, rate.cat, ...),
        future.seed = TRUE
      )

      for (i in seq_along(parallel_out)) {
        loglik[i + 1] <- parallel_out[[i]]$loglik
        model[[i + 1]] <- parallel_out[[i]]$model

      }
    
    best <- model[[which.max(loglik)]]
  })
    
  attr(best, "time") <- tt
  best
}

set.seed(427)
o1 <- simfun(nstate = 2, ntrait = 3, ntaxa = 1000, meanrate = 1)

## set multicore
multistart = 10
plan(cluster, workers = min(multistart - 1, parallel::detectCores() - 1))
oo1 <- bestfit(o1)
attr(oo1, 'time')

# ▶ Worker 1 started on PID 16380 
# ▶ Worker 2 started on PID 29912 
# ▶ Worker 3 started on PID 40796 
# ▶ Worker 4 started on PID 15692 
# ▶ Worker 5 started on PID 10716 
# ▶ Worker 6 started on PID 1040 
# ▶ Worker 7 started on PID 35284 
# ▶ Worker 8 started on PID 7112 
# ▶ Worker 9 started on PID 56600 

#   user  system elapsed 
#  26.83    0.41   60.73 

plan(sequential)
oo2 <- bestfit(o1)
attr(oo2, 'time')

# ▶ Worker 1 started on PID 22856 
# ▶ Worker 2 started on PID 22856 
# ▶ Worker 3 started on PID 22856 
# ▶ Worker 4 started on PID 22856 
# ▶ Worker 5 started on PID 22856 
# ▶ Worker 6 started on PID 22856 
# ▶ Worker 7 started on PID 22856 
# ▶ Worker 8 started on PID 22856 
# ▶ Worker 9 started on PID 22856 

#    user  system elapsed 
#  314.14    4.66  319.63 

plan(cluster, workers = min(multistart - 1, parallel::detectCores() - 1))
oo3 <- bestfit(o1, use_RTMB = TRUE)
attr(oo3, 'time')
  #  user  system elapsed 
  #  1.74    0.05   10.59 

plan(sequential)
oo4 <- bestfit(o1, use_RTMB = TRUE)
attr(oo4, 'time')
  #  user  system elapsed 
  # 15.93    0.12   26.83 