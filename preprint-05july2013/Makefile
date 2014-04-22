TARGET=emmons_drake_caps

all:
	lualatex -shell-escape $(TARGET).tex
#	bibtex $(TARGET).aux
	lualatex -shell-escape $(TARGET).tex
	sh wordcount.sh $(TARGET).tex

clean:
	rm -rf *.aux *.bbl *.blg *.log *.pdf