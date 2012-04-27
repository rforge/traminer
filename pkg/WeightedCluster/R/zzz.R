
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr<-utils:::packageDescription("WeightedCluster"))
	Tver <- TraMineR:::extract.ver(descr$Version)
	if(as.numeric(Tver[2])%%2==0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	builtDate <- strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2]
	packageStartupMessage("\nThis is ",descr$Package," ", state, " version ", descr$Version, " (Built: ", builtDate, ")")
	packageStartupMessage('\nTo get the manuals, please run:\n   vignette("WeightedCluster") ## Complete manual in French\n   vignette("WeightedClusterPreview") ## Short preview in English\n')
	packageStartupMessage("To cite WeightedCluster in publications please use:\n")
	packageStartupMessage("Studer, Matthias (2012). Étude des inégalités de genre en début de carrière\n académique à l'aide de méthodes innovatrices d'analyse de données séquentielles,\n Chapitre: Le manuel de la librairie WeightedCluster : Un guide pratique pour la\n création de typologies de trajectoires en sciences sociales avec R. Thèse de doctorat inédite, Université de Genève.\n\n")
}
