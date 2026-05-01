PDF = final.pdf
RNW = final.Rnw

.PHONY: all clean

all: $(PDF)

$(PDF): $(RNW)
	Rscript -e "knitr::knit2pdf('$(RNW)')"

clean:
	rm -f final.tex final.pdf *.aux *.log *.out *.toc *.bbl *.blg