#include "olmm.h"
#include <R_ext/Rdynload.h>
#include "Syms.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
  CALLDEF(olmm_setPar, 2),
  CALLDEF(olmm_update_marg, 2),
  CALLDEF(olmm_update_u, 1),
  CALLDEF(olmm_pred_marg, 5),
  {NULL, NULL, 0}
};

void R_init_vcrpart(DllInfo *dll)
{
  /* registering the functions */
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

  /* slots in olmm object */
  olmm_subjectSym = install("subject");
  olmm_ySym = install("y");
  olmm_XSym = install("X");
  olmm_WSym = install("W");
  olmm_weightsSym = install("weights");
  olmm_weightsSbjSym = install("weights_sbj");
  olmm_offsetSym = install("offset");
  olmm_dimsSym = install("dims");
  olmm_fixefSym = install("fixef");
  olmm_ranefCholFacSym = install("ranefCholFac");
  olmm_coefficientsSym = install("coefficients");
  olmm_etaSym = install("eta");
  olmm_uSym = install("u");
  olmm_logLikObsSym = install("logLik_obs");
  olmm_logLikSbjSym = install("logLik_sbj");
  olmm_logLikSym = install("logLik");
  olmm_scoreObsSym = install("score_obs");
  olmm_scoreSbjSym = install("score_sbj");
  olmm_scoreSym = install("score");
  olmm_infoSym = install("info");
  olmm_ghxSym = install("ghx");
  olmm_ghwSym = install("ghw");
  olmm_ranefElMatSym = install("ranefElMat");
}
