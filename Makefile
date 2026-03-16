formulas.html: formulas.qmd

%.html: %.qmd
	quarto render formulas.qmd
