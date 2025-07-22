
## our machinery for pruning algorithm via RTMB/rate matrix setup/etc.
source("../R/generate_traits.R")
source("../R/Q_template.R")
source("../R/loglik.R")

c_cache <- "corHMM_cache.rda"
if (file.exists(c_cache)) load(c_cache)
cache_save_vars <- NULL

## for 'subplex' comparison (allow ... to be passed)
remotes::install_github("astamm/nloptr", ref = remotes::github_pull("196"))
## BMB version of corHMM (for extras like returning deviance function)
remotes::install_github("bbolker/corHMM")

library(ape)
library(phangorn)
library(corHMM)
library(RTMB)
library(nloptr)

## generate a medium-size example, consistent with a
## 3-binary-trait setup

nstate <- 2; ntrait <- 3; ntaxa <- 200
set.seed(7)
m3 <- rtree(ntaxa)
g1 <- reorder(m3, "pruningwise")

Q <- setup_Q_template(n=nstate, k = ntrait)
nrates <- sum(Q!=0)
Q[Q!=0] <- rexp(nrates, rate = 10)
s <- phangorn::simSeq(g1, l = 1, Q = Q,
       type = "USER", levels = seq(nstate^ntrait),
       rate = 1)

## convert from enumerated traits (0-7) to 3 binary digits
## https://stackoverflow.com/questions/6614283/converting-decimal-to-binary-in-r
to_bin <- function(x, n = 3) {
  intToBits(x) |> rev() |> as.integer() |> tail(n)
}

setname <- function(x) {cbind(nm = rownames(x), x)}

traitMatrix <- sapply(s, function(x) to_bin(x-1)) |>
  t() |>
  as.data.frame() |>
  setname()

cache_save_vars <- c(cache_save_vars, c("t_corHMM", "fit_corHMM"))
if (!exists("fit_corHMM")) {
  t_corHMM <- system.time(
    fit_corHMM <- corHMM(phy = g1, data = traitMatrix, rate.cat = 1,
                         root.p = "maddfitz", return.devfun = TRUE)
  )
  save(list = cache_save_vars, file = c_cache)
}

## slightly redundant because we already have `s`
traitList <- multi_to_single(traitMatrix[,-1])
## trait_to_phydat <- function(x, nm) {
##   tm <- split(x, seq(x)) |> setNames(nm)
##   phangorn::phyDat(tm,type="USER", levels=unique(traitList))
## }
## dat.red <- trait_to_phydat(traitList, traitMatrix[,1])

par.score <- phangorn::parsimony(g1, s, method="fitch")/2
tl <- sum(g1$edge.length)
mean.change = par.score/tl
Q <- setup_Q_template(nstate,ntrait)
np <- sum(Q!=0)
starts <- log(sort(rexp(np, 1/mean.change), decreasing = TRUE))

Phylodata <- list(Q_template = Q,
                  tree = g1, trait_values = traitList,
                  traitMatrix = traitMatrix)
  
ff <- MakeADFun(cmb(prune_nll, Phylodata),
                list(log_trans_rates = starts), silent = TRUE)
opt.args <- NULL
t_RTMB_nlminb <- system.time(
  fit_RTMB_nlminb <- suppressWarnings(
    with(ff, do.call(nlminb, c(list(par, fn, gr), opt.args))))
)

opt.args <- list(method = "BFGS", control = list(maxit = 1000))
t_RTMB_BFGS <- system.time(
  fit_RTMB_BFGS <- suppressWarnings(
    suppressWarnings(with(ff, do.call(optim, c(list(par, fn, gr), opt.args))))
  )
)

## subplex fit as in corHMM (
opts <- list(algorithm = "NLOPT_LN_SBPLX",
             maxeval = "1000000", ftol_rel = 1.49011611938477e-08)
t_RTMB_subplex <- system.time(
  fit_RTMB_subplex <- nloptr(x0 = starts, eval_f = ff$fn, opts = opts)
)

c(corHMM = 
    -1*c(logLik(fit_corHMM)),
  RTMB_nlminb = fit_RTMB_nlminb$objective,
  RTMB_BFGS = fit_RTMB_BFGS$value,
  RTMB_subplex = fit_RTMB_subplex$objective)

## RTMB is consistently different from corHMM
##      corHMM  RTMB_nlminb    RTMB_BFGS RTMB_subplex 
##    282.5865     283.3723     283.3723     283.3772 

## try repeated starts
## mean.change from within corHMM: [1] 0.4571788
set.seed(101)
t2 <- system.time(
  r <- replicate(20, {
    ## starts <- log(sort(rexp(np, 1/mean.change), decreasing = TRUE))
    starts <- log(rexp(np, 1/mean.change))
    fit <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr)))))
    fit$objective
  })
)

plot(ecdf(r))
table(r)

## TO DO:
## use hmm2prune to compare parameters (feed corHMM parameters to RTMB)
## use devfun in BMB version of corHMM to compare (feed RTMB parameters to corHMM)

## failure to find minimum?
## difference in computing log-likelihood? (Priors, normalization constants, etc.) ?
