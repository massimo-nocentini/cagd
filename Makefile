
all:
	$(MAKE) -C bezier-deCasteljau-curves all;
	pdflatex cagd.tex

pdf:
	pdflatex cagd.tex
