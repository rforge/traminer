## .First.lib <- function(lib, pkg) {library.dynam("TraMineR", pkg, lib)}

.onAttach <- function(libname, pkgname){
	suppressWarnings(descr<-utils:::packageDescription("TraMineR"))
	packageStartupMessage("\n",descr$Package," version ", descr$Version)
	packageStartupMessage("Please visit: ", descr$URL)
	packageStartupMessage("and type 'citation(\"TraMineR\")' for information on how to cite TraMineR.\n")
}
