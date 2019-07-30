.IGNORE:

abstract = abstract
poster = bssn


abstract_pdf: TARGET=$(abstract)
abstract_pdf: pdf

abstract: TARGET=$(abstract)
abstract: ref

poster_pdf: TARGET=$(poster)
poster_pdf: pdf

poster: TARGET=$(poster)
poster: ref	

separate:
	pdftk proposal_prac.pdf cat 1-15     output main.pdf
	pdftk proposal_prac.pdf cat "16-end" output references.pdf

ref:
	pdflatex -shell-escape $(TARGET)
	bibtex    $(TARGET)
	pdflatex -shell-escape $(TARGET)
	pdflatex -shell-escape $(TARGET)
	
pdf:	
	pdflatex -shell-escape $(TARGET)


clean:
	/bin/rm -f *.dvi *.log *.aux *.toc *.bbl *blg *.end *.pdf *.ps
