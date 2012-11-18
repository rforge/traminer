

RSCRIPT="$(R_HOME)/bin/Rscript"



all: pdfclean rotateFiles WeightedCluster.pdf WeightedCluster-fr.pdf WeightedClusterPreview.pdf compactPDF clean

pdfclean:
	rm -f WeightedCluster-fr.pdf WeightedCluster.pdf WeightedClusterPreview.pdf


rotateFiles: pdfclean
	cp -f ./WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw
	cp -f ./WeightedCluster/WeightedCluster-fr-fake.Rnw WeightedCluster-fr.Rnw
	cp -f ./WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

clean: WeightedCluster.pdf WeightedCluster-fr.pdf WeightedClusterPreview.pdf
	rm -f *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary
	rm -rf graphics
	
	
WeightedClusterPreview.pdf: rotateFiles WeightedClusterPreview/WeightedClusterPreview.Rnw
	cp -f WeightedClusterPreview/WeightedClusterPreview.Rnw .
	cp -f WeightedClusterPreview/previewhctree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedClusterPreview.Rnw')"
	rm -f previewhctree.png
	cp -f WeightedClusterPreview/WeightedClusterPreview-fake.Rnw WeightedClusterPreview.Rnw

WeightedCluster.pdf: rotateFiles WeightedCluster/WeightedCluster.Rnw
	cp -f WeightedCluster/WeightedCluster.Rnw .
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	cp -f WeightedCluster/seqtree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedCluster.Rnw')"
	rm -f randB.RData wardtree.png
	cp -f WeightedCluster/WeightedCluster-fake.Rnw WeightedCluster.Rnw

WeightedCluster-fr.pdf: rotateFiles WeightedCluster/WeightedCluster-fr.Rnw
	cp -f WeightedCluster/WeightedCluster-fr.Rnw .
	cp -f WeightedCluster/randB.RData .
	cp -f WeightedCluster/wardtree.png .
	cp -f WeightedCluster/seqtree.png .
	$(RSCRIPT) --vanilla --encoding="latin1" -e "options('encoding'='latin1');library(knitr); knit2pdf('WeightedCluster-fr.Rnw')"
	rm -f randB.RData wardtree.png
	cp -f WeightedCluster/WeightedCluster-fr-fake.Rnw WeightedCluster-fr.Rnw

	
compactPDF:  WeightedCluster-fr.pdf  WeightedCluster.pdf WeightedClusterPreview.pdf
	$(RSCRIPT) --vanilla -e "library(tools); compactPDF('./', gs_quality='ebook')"