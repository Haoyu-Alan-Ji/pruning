library(Matrix)
library(ape)
library(waldo)
library(expm)
library(RTMB)

source("R/ADtools.R")
source("R/random_refit.R")

pruningAD <- function(tree, trait, state, rep.times = 100,
                      meanrate = 1, pars, opt.args = NULL) {
  Q <- Q_template(state,trait)

  repeat {
  s <- phangorn::simSeq(tree, l = 1, Q = Q,
                        type = "USER", levels = seq(nstate^ntrait),
                        rate = 1)
  if (nrow(unique(as.character(s))) == prod(nstate^ntrait)) break
  }

  s <- as.numeric(unlist(s))

  # nrates <- sum(Q != 0)
  # Q[Q!=0] <- rexp(nrates, rate = 1/meanrate)
  
  Phylodata <- list(Q_template = Q, tree = tree, trait_values = s)
  
  ff <- MakeADFun(cmb(prune_nll, Phylodata), list(log_trans_rates = pars), silent = TRUE)
  AD <- suppressWarnings(with(ff, do.call(nlminb, c(list(par, fn, gr), opt.args))))
  
  result_frame <- suppressWarnings(
    replicate(rep.times, simplify = FALSE, random_refit(ff))) |>
    do.call(what = rbind)
  
  result_good <- subset(result_frame, convergence == 0)

  if (nrow(result_good) == 0) {
    tt <- table(result_frame$convergence)
    stop("No good results found: all refits failed to converge: codes ",
         paste(sprintf("code %s = %d", names(tt), tt), collapse = ", "))
  }
  
  pars.best <- result_good[which.min(result_good$objective),
                          paste0("fitted", 1:sum(Q != 0))] |> unlist()
  gr.best <- ff$gr(pars.best)
  obj.best <- ff$fn(pars.best)
  fn.cdf <- ecdf(result_good$objective)
  
  result <- list(pars.best = pars.best,
                 obj.best = obj.best,
                 gr.best = gr.best,
                 Phylo = Phylodata,
                 pars.start = pars.start,
                 fn.cdf = fn.cdf)
  
  cat('loglik:\n'); print(obj.best)
  cat('pars:\n'); print(pars.best)
  cat('loglik cdf:\n'); plot(fn.cdf)
  class(result) <- "PAD"
  return(result)
}
