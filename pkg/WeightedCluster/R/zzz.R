
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils:::packageDescription("WeightedCluster"))
	if(utils:::packageVersion("WeightedCluster")$minor %% 2 == 0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	builtDate <- strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2]
	packageStartupMessage("This is WeightedCluster ", state, " version ", descr$Version, " (Built: ", builtDate, ")")
	packageStartupMessage('\nTo get the manuals, please run:')
	packageStartupMessage('   vignette("WeightedCluster") ## Complete manual in French')
	packageStartupMessage('   vignette("WeightedClusterPreview") ## Short preview in English')
	packageStartupMessage("\nTo cite WeightedCluster in publications please use:")
	packageStartupMessage("Studer, Matthias (2012). \u00c9tude des in\u00e9galit\u00e9s de genre en d\u00e9but de carri\u00e8re")
	packageStartupMessage("   acad\u00e9mique \u00e0 l\'aide de m\u00e9thodes innovatrices d\'analyse de donn\u00e9es")
	packageStartupMessage("   s\u00e9quentielles, Chapitre: Le manuel de la librairie WeightedCluster :")
	packageStartupMessage("   Un guide pratique pour la cr\u00e9ation de typologies de trajectoires en sciences")
	packageStartupMessage("   sociales avec R. Th\u00e8se SES 777, Facult\u00E9 des sciences \u00E9conomiques\n   et sociales, Universit\u00e9 de Gen\u00e8ve.\n")
}
