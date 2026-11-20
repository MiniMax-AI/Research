PDFLATEX = pdflatex -interaction=nonstopmode -halt-on-error

all: main.pdf

main.pdf: main.tex
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex
	$(PDFLATEX) main.tex
	rm -f main.aux main.log main.out

clean:
	rm -f main.aux main.log main.out

.PHONY: all clean
