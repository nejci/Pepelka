doc.ps:doc.dvi 
	dvips doc.dvi -o doc.ps

doc.dvi:doc.tex references.bib
	latex doc.tex && bibtex doc && latex doc.tex && latex doc.tex