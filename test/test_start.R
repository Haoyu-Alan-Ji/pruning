source("../R/generate_traits.R")
library(ape)
library(phangorn)

set.seed(7)
m3 <- rtree(20)
g1 <- reorder(m3, "pruningwise")

traitMatrix <- get_multitrait(c(state=2, trait=3), tree=g1)
traitList <- multi_to_single(traitMatrix[,-1], c(state =2 , trait = 3))
trait_to_phydat <- function(x, nm) {
  tm <- split(x, seq(x)) |> setNames(nm)
  phangorn::phyDat(tm,type="USER", levels=unique(traitList))
}
dat.red <- trait_to_phydat(traitList, traitMatrix[,1])
par.score <- phangorn::parsimony(g1, dat.red, method="fitch")/2
