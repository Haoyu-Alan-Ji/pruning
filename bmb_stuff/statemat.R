library(Matrix)

dn <- sprintf("(%d)", 1:8)
m <- matrix(c(0, 1, 2, 0, 3, 0, 0, 0, 4, 0, 0, 2, 0, 5, 0, 0, 6, 
0, 0, 1, 0, 0, 7, 0, 0, 6, 4, 0, 0, 0, 0, 8, 9, 0, 0, 0, 0, 1, 
2, 0, 0, 10, 0, 0, 4, 0, 0, 2, 0, 0, 11, 0, 6, 0, 0, 1, 0, 0, 
0, 12, 0, 6, 4, 0), nrow = 8L, ncol = 8L, dimnames = list(dn, dn))

## RTMB code:
m2 <- m
params <- (1:12)*10
m2[m2 != 0] <- params[m2[m2 != 0]]

image(Matrix(m2), colorkey = TRUE)
##
## then compute the pruning algorithm likelihood based on these rates
## likelihood of each step is based on a matrix exponential:

## x0 = {start probabilities of each state}
## M = {transition rate matrix}
## dt = {branch length/time interval for evolution}
## x1 = {end probabilities of each state} = x0 %*% expm(M * dt)

## the pruning algorithm will let you get probabilities of occupancy
## of the tips, given a starting distribution at the root, and the
## structure of the phylogeny, and the transition rate matrix

## https://github.com/thej022214/corHMM/blob/f0d1702d01ce50d86f03d9005825e6779cac9444/R/corDISC.R#L264-L276

## https://kaskr.github.io/adcomp/Introduction.html
## expm() is the matrix exponential function:
## https://kaskr.github.io/adcomp/group__matrix__functions.html#gac09966696ad9685dc8047d6c582bc0ae

## sparse matrix exponential?
## https://kaskr.github.io/adcomp/structsparse__matrix__exponential_1_1expm__series.html


## 12 parameters:
## gain/loss of SM = 2 (indep of everything else)
## gain/loss of PC = 2 (ditto)
## gain/loss of AG = {baseline, effect of PC on the rate,
##                    effect of SM on the  rate, interaction between SM/PC on AG evolution rate} x {gain/loss}

## illustration of sparse matrix: if we had a transition matrix
## like this then worrying about how to do a matrix exponential
## computation on a sparse matrix would be worth the trouble ...
d <- 100
m <- matrix(0, d, d)
m[cbind(sample(d, size = 50, replace = FALSE),
        sample(d, size = 50, replace = FALSE))] <- 1
image(Matrix(m))

## "Analyses of Phylogenetics and Evolution"
## this is the package that knows about phylogenetic tree structures
library(ape)
data("chiroptera", package = "ape")
plot(chiroptera, type = "radial")

head(chiroptera$edge)
head(reorder(chiroptera, "pruningwise")$edge)

## you should do this reordering before you start the pruning algorithm


set.seed(101)
m3 <- rtree(10)
par(mfrow=c(2,2))
head(m3$edge)
head(reorder(m3, "pruningwise")$edge)
head(reorder(m3, "postorder")$edge)
## is postorder == pruningwise in reverse??
head(reorder(m3, "postorder")$edge[18:1,])

## why do these all look identical??
plot(m3)
plot(reorder(m3, "pruningwise"))
plot(reorder(m3, "cladewise"))
plot(reorder(m3, "postorder"))

identical(m3, reorder(m3, "pruningwise"))
library(waldo)
compare(m3, reorder(m3, "pruningwise"))
