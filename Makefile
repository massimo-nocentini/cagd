
all:
	$(MAKE) -C bezier-deCasteljau-curves all;
	$(MAKE) -C b-splines all;
	pdflatex cagd.tex

pdf:
	pdflatex cagd.tex
