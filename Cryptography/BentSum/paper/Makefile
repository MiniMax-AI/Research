PDF = main.pdf
SRC = main.tex
BIB = refs.bib

.PHONY: all clean

all: $(PDF)

$(PDF): $(SRC) $(BIB)
	pdflatex -interaction=nonstopmode -halt-on-error $(SRC)
	bibtex main
	pdflatex -interaction=nonstopmode -halt-on-error $(SRC)
	pdflatex -interaction=nonstopmode -halt-on-error $(SRC)
	rm -f main.aux main.bbl main.blg main.log main.out

clean:
	rm -f main.aux main.bbl main.blg main.log main.out
