source("../R/generate_traits.R")
source("../R/Q_template.R")
source("../R/loglik.R")
library(ape)
library(phangorn)
library(corHMM)
library(RTMB)

set.seed(7)
m3 <- rtree(20)
g1 <- reorder(m3, "pruningwise")

nstate <- 2
ntrait <- 3
traitMatrix <- get_multitrait(c(nstate=2, ntrait=ntrait), tree=g1)
t1 <- system.time(corHMM_fit <- corHMM(g1, traitMatrix, rate.cat = 1,
                                       root.p = "maddfitz"))

###
traitList <- multi_to_single(traitMatrix[,-1], c(nstate = nstate ,
                                                 ntrait = ntrait))
trait_to_phydat <- function(x, nm) {
  tm <- split(x, seq(x)) |> setNames(nm)
  phangorn::phyDat(tm,type="USER", levels=unique(traitList))
}
dat.red <- trait_to_phydat(traitList, traitMatrix[,1])
par.score <- phangorn::parsimony(g1, dat.red, method="fitch")/2
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
AD <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr), opt.args))))
opt.args <- list(method = "BFGS", control = list(maxit = 1000))
AD <- suppressWarnings(with(ff, do.call(optim, c(list(par, fn, gr), opt.args))))
