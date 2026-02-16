## next thing to do: build the 'scaled posterior function' (i.e. log-prior
##  + log-likelihood) from pruning-repo machinery instead of full corHMM,
##  so that we can add in the priors ...

## copied from corhmm_logpostfun in R/utils.R

#' @param p parameters (log-hazard rates)
#' @param lb lower bound(s) for baseline priors
#' @param ub upper bound(s)
#' @param range width of Gaussian (+/- SD between mean and lower/upper bounds)
#' @param gainloss_pairs
#' @param lb_gainloss
#' @param ub_gainloss
#' @param range_gainloss number of SDs from center to lower/upper bounds
#' @param nllfun \emph{negative} log-likelihood function
#' @param negative return negative log posterior?
RTMB_logpostfun <- function(
                            ## add whatever arguments the RTMB pruning algorithm loglik function
                            ## needs (tree, trait data, etc.)
                            p,
                            lb = log(1e-9),
                            ub = log(1e2),
                            range = 3,
                            gainloss_pairs = NULL,
                            lb_gainloss = log(1e-3),
                            ub_gainloss = log(1e3),
                            range_gainloss = 3,
                            negative = FALSE
                            ) {
  prior.mean <- (lb + ub) / 2
  prior.sd <- (ub - lb) / (2 * range)
  ## call RTMB pruning-algorithm code here to compute log-likelihood ...
  loglik <- -1*nllfun(p)
  log.prior <- sum(dnorm(p, mean = prior.mean, sd = prior.sd, log = TRUE))
  ## product of likelihood and prior -> sum of LL and log-prior
  res <- loglik + log.prior
  ## calculate gain/loss priors
  if (!is.null(gainloss_pairs)) {
    gl.prior.mean <- (lb_gainloss + ub_gainloss) / 2
    gl.prior.sd <- (ub_gainloss - lb_gainloss) / (2 * range_gainloss)
    ## vapply might not work in RTMB? replace with for loop?
    gl.values <- vapply(gainloss_pairs,
                        function(x) p[x[1]] - p[x[2]],
                        FUN.VALUE = numeric(1))
    gl.log.prior <- sum(dnorm(gl.values, mean = gl.prior.mean, sd = gl.prior.sd,
                              log = TRUE))
    res <- res + gl.log.prior
  }
  if (negative) res <- -1*res
  return(res)
}


library(RTMB)
f <- function(x) x^2
g <- function(p) {
    getAll(p)
    f(x)^2
}

ff <- MakeADFun(g, parameters = list(x=1))
ff$fn(1)
ff$fn(2)
