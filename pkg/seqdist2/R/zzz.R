# Pierre-Alexandre Fonta (2016-2017)

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(strwrap("<!> seqdist2 is now deprecated in favour of
    TraMineR 2.x functions seqdist() and seqcost()"))
}
