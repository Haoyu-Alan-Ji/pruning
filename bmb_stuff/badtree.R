library(ape)
edge1 <- matrix(c(
  5, 3,
  3, 1,
  3, 2,
  5, 4
), byrow = TRUE, ncol = 2)

tip.label <- c("t1", "t2", "t3")

edge.length <- rep(1, nrow(edge1))

Nnode <- 2

tree <- list(edge = edge1,
             tip.label = tip.label,
             edge.length = edge.length,
             Nnode = Nnode)
class(tree) <- "phylo"

## plot(tree, type = "phylogram", show.node.label = TRUE,  edge.width = 2, cex = 0.8)
## tiplabels()
## nodelabels()
tree$edge <- matrix(c(
  5, 1,
  5, 4,
  4, 2,
  4, 3),
  byrow = TRUE, ncol = 2)
## plot(tree)  ## core dump!

str(tree)
is.binary(tree)
is.rooted(tree)
is.ultrametric(tree)

set.seed(101)
tree2 <- rtree(3)
unclass(tree2)

## edge: a two-column matrix of mode numeric where each row represents
##       an edge of the tree; the nodes and the tips are symbolized
##       with numbers; the tips are numbered 1, 2, ..., and the nodes
##       are numbered after the tips. For each row, the first column
##       gives the ancestor.

plot(tree2)
tiplabels()
nodelabels()

