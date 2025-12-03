##' @examples
#' realfun(raw_traitM = c_traits, raw_tree = c_tree)
realfun <- function(nstate = 2, raw_traitM, raw_tree) {
  trait_number <- ncol(raw_traitM) - 1L
  ss <- function() {
    sample(0:(nstate - 1L), size = Ntip(raw_tree), replace = TRUE)
  }
  repeat {
    d <- as.data.frame(replicate(trait_number, ss(), simplify = FALSE))
    names(d) <- paste0("trait", seq_len(trait_number))
    if (nrow(unique(d)) == nstate^trait_number) break
  }
  d <- cbind(Species = raw_tree$tip.label, d)
  rownames(d) <- NULL
  list(tree = raw_tree, data = d)
}

##' @examples
#' o1 <- sample_subtree(c_tree, c_traits, 50)
#' corHMM(phy = o1$tree, data = o1$traits, rate.cat = 1, use_RTMB = TRUE)
sample_subtree <- function(tree, traits, ntaxa) {
  
  if (ntaxa > Ntip(tree)) {
    stop("ntaxa is greater than nodes")
  }
  
  keep_tips <- sample(tree$tip.label, ntaxa)
  
  sub_tree <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))
  
  sub_traits <- traits[traits$Species %in% sub_tree$tip.label, ]
  sub_traits <- sub_traits[match(sub_tree$tip.label, sub_traits$Species), ]
  
  stopifnot(all(sub_tree$tip.label == sub_traits$Species))
  
  list(tree = sub_tree, traits = sub_traits)
}