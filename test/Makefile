all: p1_simtrait.html g2_sim_8by8.html

%.html: %.qmd
	quarto render $<

p1_simtrait.html: p1_simtrait.qmd ../R/pruning_funs.R
	quarto render $<

g1_sim_4by4.html: g1_sim_4by4.qmd ../R/pruning_funs.R
	quarto render $<

g2_sim_8by8.html: g2_sim_8by8.qmd ../R/pruning_funs.R
	quarto render $<
