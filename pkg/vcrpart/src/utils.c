#include "utils.h"
#include <Rmath.h>

/**
 * ---------------------------------------------------------
 * Duplicate an R object
 *
 * @param x an R object.
 *
 * @return the duplicate of 'x'
 * ---------------------------------------------------------
 */

SEXP vcrpart_duplicate(SEXP x) {
  SEXP rval = PROTECT(duplicate(x));
  UNPROTECT(1);
  return rval;
}

/**
 * ---------------------------------------------------------
 * Computes the split matrix for nominal variables for
 * the 'tvcm' alorithm.
 *
 * @param xlev integer vector, all observed categories.
 * @param nl total number of categories.
 * @param xdlev integer vector, the levels observed.
 *    in the current node.
 * @param nld number of categories observed in the current
 *    node.
 * @param mi number of possible splits with the categories
 *    observed in the current node.
 *
 * @return an integer vector of length nld * mi to be used
 *    to construct the split matrix.
 * ---------------------------------------------------------
 */

SEXP tvcm_nomsplits(SEXP nl, SEXP xdlev, SEXP nld, SEXP mi) {
  
  int *xdlev_c = INTEGER(xdlev);
  const int nl_c = INTEGER(nl)[0], 
    nld_c = INTEGER(nld)[0], 
    mi_c = INTEGER(mi)[0];

  SEXP indx = PROTECT(allocVector(INTSXP, nl_c * mi_c));
  for (int i = 0; i < (nl_c * mi_c); i++) INTEGER(indx)[i] = 0;
  int ii = 0;
  if (mi_c > 0) {
    for (int i = 0; i < mi_c; i++) {
      ii = i + 1;
      for (int l = 0; l < nld_c; l++) {
	INTEGER(indx)[i * nl_c + xdlev_c[l] - 1] = ii % 2;
	ii = ii / 2;
      }
    }
  }
  UNPROTECT(1);
  return indx;
}
