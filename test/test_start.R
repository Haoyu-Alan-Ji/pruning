## our machinery for pruning algorithm via RTMB/rate matrix setup/etc.
source("../R/generate_traits.R")
source("../R/Q_template.R")
source("../R/loglik.R")

c_cache <- "corHMM_cache.rda"
if (file.exists(c_cache)) load(c_cache)
cache_save_vars <- NULL

## for 'subplex' comparison (allow ... to be passed)
try(remotes::install_github("astamm/nloptr", ref = remotes::github_pull("196")))
## BMB version of corHMM (for extras like returning deviance function)
try(remotes::install_github("bbolker/corHMM"))

library(ape)
library(phangorn)
library(corHMM)
library(RTMB)
library(nloptr)
library(Matrix)
library(gridExtra)
library(ggplot2); theme_set(theme_bw())  ## try tinyplot?

## generate a medium-size example, consistent with a
## 3-binary-trait setup

nstate <- 2; ntrait <- 3; ntaxa <- 200
set.seed(7)
m3 <- rtree(ntaxa)
g1 <- reorder(m3, "pruningwise")

Q0 <- Q <- setup_Q_template(n=nstate, k = ntrait)
nrates <- sum(Q!=0)
Q[Q!=0] <- rexp(nrates, rate = 10)
## simulate traits
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

## corHMM fit
cache_save_vars <- c(cache_save_vars, c("t_corHMM", "fit_corHMM"))
if (!exists("fit_corHMM")) {
  t_corHMM <- system.time(
    fit_corHMM <- corHMM(phy = g1, data = traitMatrix, rate.cat = 1,
                         root.p = "maddfitz", return.devfun = TRUE)
  )
  save(list = cache_save_vars, file = c_cache)
}

## trait coding differs?

## 4 - 7
## 2 - 5
## 6 - 6
## 1 - 1
## 5 - 2
## 3 - 3
## 8 - 8
perm <- c(1,5,3,7,2,6,4,8)
## so (corHMM states)[perm] → (RTMB states)
##  RTMB states[order(perm)] -> (corHMM states)
Qperm <- Q0[order(perm), order(perm)]
rperm <- Qperm[Qperm!=0]

chk <- function() print(rperm)
chk()


## slightly redundant because we already have `s`
traitList <- multi_to_single(traitMatrix[,-1])

## if we wanted to convert *from* traits to phangorn-style phyDat
## trait_to_phydat <- function(x, nm) {
##   tm <- split(x, seq(x)) |> setNames(nm)
##   phangorn::phyDat(tm,type="USER", levels=unique(traitList))
## }
## dat.red <- trait_to_phydat(traitList, traitMatrix[,1])

## compute parsimony score, to set scale for initial conditions
par.score <- phangorn::parsimony(g1, s, method="fitch")/2
tl <- sum(g1$edge.length)
mean.change = par.score/tl
Q <- setup_Q_template(nstate,ntrait)
np <- sum(Q!=0)
starts <- log(sort(rexp(np, 1/mean.change), decreasing = TRUE))

chk()

## set up RTMB objective function
Phylodata <- list(Q_template = Q,
                  tree = g1, trait_values = traitList,
                  traitMatrix = traitMatrix)
ff <- MakeADFun(cmb(prune_nll, Phylodata),
                list(log_trans_rates = starts), silent = TRUE)

chk()

## various fits ...
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

## subplex fit as in corHMM
opts <- list(algorithm = "NLOPT_LN_SBPLX",
             maxeval = "1000000", ftol_rel = 1.49011611938477e-08)
t_RTMB_subplex <- system.time(
  fit_RTMB_subplex <- nloptr(x0 = starts, eval_f = ff$fn, opts = opts)
)

chk()

c(corHMM = 
    -1*c(logLik(fit_corHMM)),
  RTMB_nlminb = fit_RTMB_nlminb$objective,
  RTMB_BFGS = fit_RTMB_BFGS$value,
  RTMB_subplex = fit_RTMB_subplex$objective)

## RTMB is consistent across optimizers, but different from corHMM
##      corHMM  RTMB_nlminb    RTMB_BFGS RTMB_subplex 
##    282.5927     283.3723     283.3723     283.3772


cbind(t_corHMM,
      t_RTMB_nlminb,
      t_RTMB_BFGS,
      t_RTMB_subplex)[c("user.self","sys.self","elapsed"),]

##           t_corHMM t_RTMB_nlminb t_RTMB_BFGS t_RTMB_subplex
## user.self  171.714          0.98       0.686         10.593
## sys.self   291.099          0.00       0.000          0.000
## elapsed    116.101          0.98       0.687         10.596

chk()

## try repeated starts
## mean.change from within corHMM: [1] 0.4571788
set.seed(101)
t2 <- system.time(
  r <- replicate(20, {
    cat(".")
    ## starts <- log(sort(rexp(np, 1/mean.change), decreasing = TRUE))
    starts <- log(rexp(np, 1/mean.change))
    fit <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr)))))
    fit$objective
  })
)
## all identical, regardless of starting value
table(r)
plot(ecdf(r))
table(r)

chk()

## index matrix from fit_corHMM
## non-zero/NA positions are the same
stopifnot(identical(which(!is.na(fit_corHMM$index.mat), arr.ind = TRUE),
                    which(Q != 0, arr.ind = TRUE)))
## order/sequence of rates is the same
stopifnot(identical(c(na.omit(c(fit_corHMM$index.mat))), Q[Q!=0]))

chk()

Q_sol <- Q
Q_sol[Q_sol!=0] <- exp(fit_RTMB_nlminb$par)

Q_sol_corHMM <- fit_corHMM$solution
Q_sol_corHMM[is.na(Q_sol_corHMM)] <- 0

m1 <- Matrix::Matrix(Q_sol)
m2 <- Matrix::Matrix(Q_sol_corHMM)[perm, perm]
print(m1, digits = 3)
print(m2, digits = 3)

grid.arrange(image(m1), image(m2), nrow = 1)

chk()

## compare computed log-likelihoods

print(rperm)  ## how did rperm get corrupted?

with(fit_corHMM, do.call(devfun, args.list)) ## 282.5927, recover value from fit
ff$fn(fit_corHMM$args.list$p[order(rperm)]) ## 283.4434

alist <- fit_corHMM$args.list
alist$p <- fit_RTMB_nlminb$par[rperm]
do.call(fit_corHMM$devfun, alist)  ## 282.6504
ff$fn(fit_RTMB_nlminb$par) ## 283.3723, recover value from fit

all.equal(unname(m1), unname(m2))  ## 5% mean relative difference (on rate scale, not log-rate scale)

## conclusion: log-likelihood calcs differ slightly (numerical fuzz?)
## corHMM prefers its estimates over RTMB's by 0.06 log-likelihood units
## RTMB prefers its estimates over corHMM's by 0.07 log-likelihood units

## is there a good metric for how much this matters in some more concrete terms?
## would like to use RTMB::sdreport(ff) but (probably) the
H <- optimHess(fit_RTMB_nlminb$par, ff$fn, ff$gr)
## add ridge penalty so we can invert ...
V <- solve(H + 1e-5*diag(nrow(H)))
RTMB_sd <- sqrt(diag(V))
RTMB_sd[RTMB_sd>300] <- NA
nrates <- length(RTMB_sd)
plot_dat <- data.frame(method = rep(c("RTMB", "corHMM"), each = nrates),
                       par = factor(rep(1:nrates, 2)),
                       est = c(fit_RTMB_nlminb$par[rperm], fit_corHMM$args.list$p),
                       se = c(RTMB_sd, rep(NA, nrates)))
plot_dat <- transform(plot_dat, par = reorder(par, est))
ggplot(plot_dat, aes(est, par, colour = method)) +
  geom_pointrange(aes(xmin = est-2*se, xmax = est+2*se), position = position_dodge(width = 0.5))








