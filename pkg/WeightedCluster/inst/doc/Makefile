

RSCRIPT="$(R_HOME)/bin/Rscript"


all: pdfclean rotateFiles WeightedClusterPreview.pdf clean

pdfclean:
	rm -f WeightedClusterPreview.pdf


rotateFiles:
##	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw
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
