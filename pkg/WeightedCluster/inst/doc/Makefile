

RSCRIPT="$(R_HOME)/bin/Rscript"



all: pdfclean rotateFiles WeightedCluster.pdf WeightedCluster-fr.pdf WeightedClusterPreview.pdf clean

pdfclean:
	rm -f WeightedCluster-fr.pdf WeightedCluster.pdf WeightedClusterPreview.pdf


rotateFiles:
	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw
	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster-fr.Rnw
	cp -f ./WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

clean:
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary
	rm -rf graphics
	
	
WeightedClusterPreview.pdf: WeightedClusterPreview/WeightedClusterPreview.Rnw
	cp -f WeightedClusterPreview/WeightedClusterPreview.Rnw .
	cp -f WeightedClusterPreview/previewhctree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedClusterPreview.Rnw')"
	rm -f previewhctree.png
	cp -f WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

WeightedCluster.pdf: WeightedCluster/WeightedCluster.Rnw
	cp -f WeightedCluster/WeightedCluster.Rnw .
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedCluster.Rnw')"
	rm -f randB.RData wardtree.png
	cp -f WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw

WeightedCluster-fr.pdf: WeightedCluster/WeightedCluster-fr.Rnw
	cp -f WeightedCluster/WeightedCluster-fr.Rnw .
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedCluster-fr.Rnw')"
	rm -f randB.RData wardtree.png
	cp -f WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster-fr.Rnw

	
compactPDF:  WeightedCluster-fr.pdf  WeightedCluster.pdf WeightedClusterPreview.pdf
	$(RSCRIPT) --vanilla -e "library(tools); compactPDF('./', gs_quality='ebook')"