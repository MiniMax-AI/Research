MAIN = main

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(MAIN).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(MAIN).tex
	rm -f $(MAIN).aux $(MAIN).log $(MAIN).out
clean:
	rm -f $(MAIN).aux $(MAIN).log $(MAIN).out $(MAIN).toc $(MAIN).synctex.gz

.PHONY: all clean
