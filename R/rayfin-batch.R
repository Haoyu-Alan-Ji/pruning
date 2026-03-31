
library(here)
## ----pkgs, include=FALSE, message = FALSE, warning=FALSE-------------------------------------------------------------------------------------
library(Matrix)
library(ape)
library(tidyverse); theme_set(theme_bw())
library(expm)
library(RTMB)
library(colorspace)


## ----seed, include=FALSE---------------------------------------------------------------------------------------------------------------------
set.seed(427)


## ----funs, message = FALSE, warning=FALSE----------------------------------------------------------------------------------------------------
#| code-fold: false

## corHMM version
## https://github.com/bbolker/corHMM/blob/bolker_clean/misc/test_RTMB.R
## remotes::install_github("bbolker/corHMM", ref = "bolker_clean")
source(here::here("R", "pruning_funs.R"))
source(here::here("R", "benchmark.R"))
source(here::here("R", "realtree.R"))


## ----realdat---------------------------------------------------------------------------------------------------------------------------------
#| code-summary: "import data"
fish_traits <- read.csv(here::here('data', "binaryTraitData.csv"))
fish_tree   <- readRDS(here::here('data', "treeSingle.rds"))


## ----clean3traits----------------------------------------------------------------------------------------------------------------------------
#| code-summary: "data cleaning"
ctraits <- fish_traits[fish_traits$species %in% fish_tree$tip.label, ]

ctraits <- ctraits[match(fish_tree$tip.label, ctraits$species), ]

c_traits <- data.frame(
  Species = ctraits$species,
  ag = ctraits$ag,
  care = ctraits$care,
  spawning = ctraits$spawning
)

c_traits <- c_traits[!is.na(c_traits$care) & !is.na(c_traits$spawning), ]
c_tree <- drop.tip(fish_tree, setdiff(fish_tree$tip.label, c_traits$Species))

c_traits <- c_traits[match(c_tree$tip.label, c_traits$Species), ]


## ----getdat----------------------------------------------------------------------------------------------------------------------------------
#| code-summary: "resampling and running"
ntaxvec <- as.numeric(45:234)

dd <- expand.grid(seed = 101:104, ntaxa  = ntaxvec, ntrait = 3, model = "ARD", rate.cat = 1) |>
  transform(model = as.character(model))

res <- list()

for (i in seq_len(nrow(dd))) {
  cat(sprintf("\rRunning %d/%d", i, nrow(dd)))

  ntaxa_i <- dd$ntaxa[i]
  sub <- sample_subtree(tree = c_tree, traits = c_traits, ntaxa  = ntaxa_i)

  args_i <- c(
    as.list(dd[i, ]),   
    list(
      traitMatrix = sub$traits, 
      realtree = sub$tree 
    )
  )

  res[[i]] <- do.call(sumfun, args_i)
}

res_df <- do.call(rbind, res)
saveRDS(res_df, file = here::here("data", "trait3_subtree.rds"))

