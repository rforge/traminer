

RSCRIPT="$(R_HOME)/bin/Rscript"


all: pdfclean rotateFiles WeightedCluster.pdf WeightedClusterPreview.pdf clean

pdfclean:
	rm -f WeightedCluster.pdf WeightedClusterPreview.pdf


rotateFiles:
	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw
	cp -f ./WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary
	rm -rf graphics
	
	
WeightedClusterPreview.pdf: WeightedClusterPreview/WeightedClusterPreview.Rnw
	cp -f WeightedClusterPreview/WeightedClusterPreview.Rnw .
	cp -f WeightedClusterPreview/previewhctree.png .
	$(RSCRIPT) --vanilla -e "library(knitr); knit2pdf('WeightedClusterPreview.Rnw');purl('WeightedClusterPreview.Rnw', output='WeightedClusterPreview.R')"
	rm -f previewhctree.png
	cp -f WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

WeightedCluster.pdf: WeightedCluster/WeightedCluster.Rnw
	cp -f WeightedCluster/WeightedCluster.Rnw .
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	$(RSCRIPT) --vanilla -e "library(knitr); knit2pdf('WeightedCluster.Rnw')"
	## $(RSCRIPT) --vanilla -e "tools::texi2dvi( 'WeightedCluster.tex', pdf = TRUE, clean = FALSE )"
	## bibtex WeightedCluster
	## $(RSCRIPT) -e "tools::texi2dvi( 'WeightedCluster.tex', pdf = TRUE, clean = TRUE )"
	rm -f randB.RData wardtree.png
	cp -f WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw

compactPDF:  WeightedCluster.pdf WeightedClusterPreview.pdf
	$(RSCRIPT) --vanilla -e "library(tools); compactPDF('./', gs_quality='ebook')"