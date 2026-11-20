PDF = main.pdf
SRC = main.tex

.PHONY: all clean

all: $(PDF)

$(PDF): $(SRC)
	pdflatex -interaction=nonstopmode -halt-on-error $(SRC)
	pdflatex -interaction=nonstopmode -halt-on-error $(SRC)
	rm -f main.aux main.log main.out

clean:
	rm -f main.aux main.log main.out
