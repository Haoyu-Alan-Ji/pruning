##' @example 
##' if (FALSE) {
##' system.time(Q_big <- setup_Q_template(k=5)) ## 10 seconds
##' png("bigmat.png")
##' imat(Q_big)
##' dev.off()
##' }
imat <- function(m, useRaster = TRUE, ...) {
  require(Matrix)
  image(Matrix(m), useRaster = useRaster, xlab = "", ylab = "", sub = "", ...)
}

## thought I could do a Kronecker product trick but now I don't see how ...
## if (use_kron) {
##   Q0 <- function(n) {
##     m <- matrix(1, n, n)
##     diag(m) <- 0
##     m
##   }
##   return(Reduce(kronecker, lapply(n, Q0)))
## }