#include "olmm.h"
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

extern
#include "Syms.h"

#ifdef __GNUC__
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
/* this is necessary (and sufficient) for Solaris 10: */
# ifdef __sun
#  include <alloca.h>
# endif
#endif

/* allocate n elements of type t */
#define Alloca(n, t)  (t *) alloca( (size_t) ( (n) * sizeof(t) ) )

/* set array values to zero */
#define AllocVal(x, n, v) {int _I_, _SZ_ = (n); for(_I_ = 0; _I_ < _SZ_; _I_++) (x)[_I_] = (v);}

/* positions in the dims vector */
enum dimP {
  n_POS = 0,      /* number of observations */
    N_POS,        /* number of subject clusters */
    p_POS,        /* number of fixed effects */ 
    pEta_POS,     /* number of columns of X */
    pInt_POS,     /* number of colums that takes the intercept in X */
    pEtaVar_POS,  /* number of predictor-variable fixef */
    pEtaInv_POS,  /* number of predictor-invariant fixef */
    q_POS,        /* number of random coefficients */
    qEta_POS,     /* number of columns of W */
    qEtaVar_POS,  /* number of predictor-variable ranef */
    qEtaInv_POS,  /* number of predictor-invariant ranef */
    J_POS,        /* number of ordinal response levels */
    nEta_POS,     /* number of predictors = J - 1 */
    nPar_POS,     /* total number of parameters */
    nGHQ_POS,     /* number of quadrature points per dimension */
    nQP_POS,      /* number of quadrature points total */
    fTyp_POS,     /* model family (1 = cumulative, 2 = ...) */
    lTyp_POS,     /* link function (1 = logit, 2 = ...)*/
    verb_POS,     /* verbose */
    numGrad_POS,  /* whether scores are computed numerically */
    numHess_POS   /* whether Hessian is computed numerically */
    };

static R_INLINE double *SLOT_REAL_NULL(SEXP obj, SEXP nm)
{
  SEXP pt = GET_SLOT(obj, nm);
  return LENGTH(pt) ? REAL(pt) : (double*) NULL;
}

#define Y_SLOT(x) INTEGER(GET_SLOT(x, olmm_ySym))
#define X_SLOT(x) SLOT_REAL_NULL(x, olmm_XSym)
#define W_SLOT(x) SLOT_REAL_NULL(x, olmm_WSym)
#define WEIGHTS_SLOT(x) SLOT_REAL_NULL(x, olmm_weightsSym)
#define WEIGHTSSBJ_SLOT(x) SLOT_REAL_NULL(x, olmm_weightsSbjSym)
#define OFFSET_SLOT(x) SLOT_REAL_NULL(x, olmm_offsetSym)
#define DIMS_SLOT(x) INTEGER(GET_SLOT(x, olmm_dimsSym))
#define FIXEF_SLOT(x) SLOT_REAL_NULL(x, olmm_fixefSym)
#define RANEFCHOLFAC_SLOT(x) SLOT_REAL_NULL(x, olmm_ranefCholFacSym)
#define COEFFICIENTS_SLOT(x) SLOT_REAL_NULL(x, olmm_coefficientsSym)
#define ETA_SLOT(x) SLOT_REAL_NULL(x, olmm_etaSym)
#define U_SLOT(x) SLOT_REAL_NULL(x, olmm_uSym)
#define LOGLIKOBS_SLOT(x) SLOT_REAL_NULL(x, olmm_logLikObsSym)
#define LOGLIKSBJ_SLOT(x) SLOT_REAL_NULL(x, olmm_logLikSbjSym)
#define LOGLIK_SLOT(x) SLOT_REAL_NULL(x, olmm_logLikSym)
#define SCOREOBS_SLOT(x) SLOT_REAL_NULL(x, olmm_scoreObsSym)
#define SCORESBJ_SLOT(x) SLOT_REAL_NULL(x, olmm_scoreSbjSym)
#define SCORE_SLOT(x) SLOT_REAL_NULL(x, olmm_scoreSym)
#define INFO_SLOT(x) SLOT_REAL_NULL(x, olmm_infoSym)
#define GHX_SLOT(x) SLOT_REAL_NULL(x, olmm_ghxSym)
#define GHW_SLOT(x) SLOT_REAL_NULL(x, olmm_ghwSym)
#define RANEFELMAT_SLOT(x) SLOT_REAL_NULL(x, olmm_ranefElMatSym)

/**
 * ---------------------------------------------------------
 * Utility functions (may relocate them ...)
 * ---------------------------------------------------------
 */

/* set link function */

double olmm_gLink(double x, int link) {
  switch (link) {
  case 1: // logit
    return dlogis(x, 0.0, 1.0, 0.0);
  case 2: // probit
    return dnorm(x, 0.0, 1.0, 0.0);
  case 3: // cauchit
    return dcauchy(x, 0.0, 1.0, 0.0);
  default : // all other
    error("link not recognised\n");
    return NA_REAL;
  }
}

double olmm_GLink(double x, int link) {
  switch (link) {
  case 1: // logit
    return plogis(x, 0.0, 1.0, 1.0, 0.0);
  case 2: // probit
    return pnorm(x, 0.0, 1.0, 1.0, 0.0);
  case 3: // cauchit
    return pcauchy(x, 0.0, 1.0, 1.0, 0.0);
  default : // all other
    error("link not recognised\n");
    return NA_REAL;
  }
}


/**
 * ---------------------------------------------------------
 * Store a new parameter vector in a olmm object
 *
 * @param x a olmm object
 * @param par the new parameter vector to store
 * ---------------------------------------------------------
 */

SEXP olmm_setPar(SEXP x, SEXP par) {
  
  int *dims = DIMS_SLOT(x);

  double *fixef = FIXEF_SLOT(x), /* used slots */
    *vecRanefCholFac = RANEFCHOLFAC_SLOT(x),
    *coefficients = COEFFICIENTS_SLOT(x), 
    *ranefElMat = RANEFELMAT_SLOT(x),
    *newPar = REAL(par);

  const int p = dims[p_POS], family = dims[fTyp_POS],  
    pEtaVar = dims[pEtaVar_POS], pEtaInv = dims[pEtaInv_POS],
    pEta = dims[pEta_POS], q = dims[q_POS], J = dims[J_POS], 
    nEta = dims[nEta_POS], nPar = dims[nPar_POS],
    lenVecRanefCholFac = q * q, lenVRanefCholFac = q * (q + 1) / 2,
    i1 = 1;

  const double zero = 0.0, one = 1.0; /* for matrix manipulations */
  
  double *tmpVec = Alloca(pEtaInv, double);

  /* overwrite coefficients slot */
  Memcpy(coefficients, newPar, nPar);

  /* overwrite fixed effect parameter slot */
  for (int i = 0; i < nEta; i++) {
    Memcpy(fixef + i * pEta, newPar + i * pEtaVar, pEtaVar);
    Memcpy(fixef + i * pEta + pEtaVar, newPar + pEtaVar*nEta, pEtaInv);
  }

  /* set predictor-invariant fixed effects of adjacent-category model 
     to use the Likelihood of the baseline-category model */
  if (family == 3) {
    for (int i = 0; i < pEtaInv; i++) {
      for (int j = 0; j < nEta; j++) {
	fixef[(pEtaVar + pEtaInv) * j + pEtaVar + i] *=
	  J - j - 1;
      }
    }
  }

  /* overwrite (cholesky decompositioned) random effect parameters */
  double *vRanefCholFac = Memcpy(Alloca(lenVRanefCholFac, double), 
				 newPar + p, lenVRanefCholFac);

  F77_CALL(dgemv)("T", /* multiply with l. t. elimination matrix */ 
		  &lenVRanefCholFac, &lenVecRanefCholFac,
  		  &one, ranefElMat, 
		  &lenVRanefCholFac, vRanefCholFac, &i1, 
		  &zero, vecRanefCholFac, &i1);

  return R_NilValue;
  }


/**
 * ---------------------------------------------------------
 * Update marginal log-Likelihood, score and info slots.
 *    
 * @param  x a olmm object
 *
 * @return R_NilValue
 * ---------------------------------------------------------
 */

SEXP olmm_update_marg(SEXP x, SEXP par) {

  /* get subject slot */  
  SEXP subject_0 = GET_SLOT(x, olmm_subjectSym), subject_1;
  PROTECT(subject_1 = coerceVector(subject_0, INTSXP));

  /* get targed variable slot */
  SEXP y_0 = GET_SLOT(x, olmm_ySym), y_1;
  PROTECT(y_1 = coerceVector(y_0, INTSXP));

  /* integer valued slots and pointer to factor valued slots */
  int *y = INTEGER(y_1), *subject = INTEGER(subject_1),
    *dims = DIMS_SLOT(x);
  R_CheckStack();

  /* numeric valued objects */
  double *X = X_SLOT(x), *W = W_SLOT(x), 
    *weights_obs = WEIGHTS_SLOT(x), *weights_sbj = WEIGHTSSBJ_SLOT(x),
    *offset = OFFSET_SLOT(x), *eta = ETA_SLOT(x), 
    *fixef = FIXEF_SLOT(x), *ranefCholFac = RANEFCHOLFAC_SLOT(x),
    *logLik_sbj = LOGLIKSBJ_SLOT(x), *logLik = LOGLIK_SLOT(x), 
    *score_obs = SCOREOBS_SLOT(x), *score_sbj = SCORESBJ_SLOT(x), 
    *score = SCORE_SLOT(x),
    *info = INFO_SLOT(x),
    *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x), 
    *ranefElMat = RANEFELMAT_SLOT(x);
  R_CheckStack();

  /* set constants (dimensions of vectors etc.) */
  const int n = dims[n_POS], N = dims[N_POS], 
    p = dims[p_POS], pEta = dims[pEta_POS],
    pEtaVar = dims[pEtaVar_POS], pEtaInv = dims[pEtaInv_POS], 
    q = dims[q_POS], qEta = dims[qEta_POS], 
    qEtaVar = dims[qEtaVar_POS], qEtaInv = dims[qEtaInv_POS],
    J = dims[J_POS], nPar = dims[nPar_POS], 
    nEta = dims[nEta_POS], nGHQ = dims[nGHQ_POS], nQP = dims[nQP_POS], 
    family = dims[fTyp_POS], link = dims[lTyp_POS], 
    verb = dims[verb_POS], numGrad = dims[numGrad_POS],
    numHess = dims[numHess_POS],
    lenVecRanefCholFac = q * q, lenVRanefCholFac = q * (q + 1) / 2;

  /* variables for matrix operations etc. */
  int i1 = 1, tmpJ, subsTmp; 
  double one = 1.0, zero = 0.0, 
    gq_weight = 1, scoreCondInv, scTermBL, 
    logLikCond_modified;
      
  /* define internal objects */
  double *etaCLM = (double*) NULL,
    *etaRanefCLM = (double*) NULL,
    *sumBL = (double*) NULL,
    *etaRanef = Calloc(n * nEta, double),
    *ranefTmp = Alloca(qEtaInv, double),
    *gq_nodes = Alloca(q, double),
    *ranefVec = Alloca(q, double),
    *ranef = Alloca(qEta * nEta, double),
    *logLikCond_obs = Calloc(n, double),
    *logLikCond_sbj = Calloc(N, double),
    *scoreCondVar = Alloca(nEta, double),
    *vecRanefTerm = Alloca(lenVecRanefCholFac, double),
    *vRanefTerm = Alloca(lenVRanefCholFac, double),
    *etaTmp = Alloca(nEta, double),
    *scoreCond_obs = (double*) NULL,
    *scoreCond_sbj = (double*) NULL;
  R_CheckStack();

  if (numGrad == 0) {
    scoreCond_obs = Calloc((1 - numGrad) * n * nPar + numGrad, double);
    scoreCond_sbj = Calloc((1 - numGrad) * N * nPar + numGrad, double);
  }

  AllocVal(etaRanef, n * nEta, zero);
  
  /* update parameters */
  olmm_setPar(x, par);
  
  /* initialize eta = offset */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < nEta; j++) eta[n * j + i] = offset[i];
  
  /* update eta = offset + Xb */
  F77_CALL(dgemm)("N", "N", &n, &nEta, &pEta, &one, X, &n, fixef, 
		  &pEta, &one, eta, &n);

  switch (family) {
  case 1: 
    etaCLM = Calloc(n * 2 , double);
    etaRanefCLM = Calloc(n * 2 , double);
    for (int i = 0; i < n; i++) {
      etaCLM[i] /* lower */
	= y[i] > 1 ? eta[n * (y[i] - 2) + i] : -DBL_MAX; 
      etaCLM[i + n] = /* upper */
	y[i] < J ? eta[n * (y[i] - 1) + i] : DBL_MAX;
    }
    break;
  case 2: case 3:
    sumBL = Calloc(n, double);
    break;
  }
  R_CheckStack();

  /* initializations */
  logLik[0] = 0;
  AllocVal(logLik_sbj, N, zero);

  if (numGrad == 0) {

    AllocVal(score, nPar, zero);
    AllocVal(score_sbj, N * nPar, zero);
    AllocVal(score_obs, n * nPar, zero);
  }

  /* Gauss-Hermite quadrature */

  for (int k = 0; k < nQP; k++) {
    
    /* clear temporary objects */
    AllocVal(logLikCond_obs, n, zero);
    AllocVal(logLikCond_sbj, N, zero);

    AllocVal(scoreCond_obs, (1 - numGrad) * n * nPar + numGrad, zero);
    AllocVal(scoreCond_sbj, (1 - numGrad) * N * nPar + numGrad, zero);

    if (family == 2 | family == 3) AllocVal(sumBL, n, zero);

    /* prepare integration nodes and weights */    
    for (int i = 0; i < q; i++) {
      gq_nodes[i] = ghx[nQP * i + k];
      gq_weight *= ghw[nQP * i + k];
    }
    gq_weight = log(gq_weight);

    /* multiply ranefCholFac with actual nodes */
    F77_CALL(dgemv)("N", &q, &q, &one, ranefCholFac, &q, 
		    gq_nodes, &i1, &zero, ranefVec, &i1);

    /* ranefVec (vector) to ranef (matrix) */
    for (int i = 0; i < nEta; i++) {
      Memcpy(ranef + i * qEta, ranefVec + i * qEtaVar, qEtaVar);
      Memcpy(ranef + i*qEta + qEtaVar, ranefVec + qEtaVar*nEta, qEtaInv);
    }

    /* set predictor-invariant random effects of adjacent-category models 
       to use the Likelihood of the baseline-category model */
    if (family == 3) {
      for (int i = 0; i < qEtaInv; i++) {
	for (int j = 0; j < nEta; j++) {
	  ranef[(qEtaVar + qEtaInv) * j + qEtaVar + i] *=
	    J - j - 1; 
	}
      }
    }

    /* compute contribution of random effects to linear predictor */
    F77_CALL(dgemm)("N", "N", &n, &nEta, &qEta, &one, W, &n, 
		    ranef, &qEta, &zero, etaRanef, &n);    

    /* use etaRanefCLM to accelerate computations */
    if (family == 1) {
      for (int i = 0; i < n; i++) {
	etaRanefCLM[i] /* lower */
	  = y[i] > 1 ? etaRanef[n * (y[i] - 2) + i] : -DBL_MAX; 
	etaRanefCLM[n + i] = /* upper */
	  y[i] < J ? etaRanef[n * (y[i] - 1) + i] : DBL_MAX;
      }
    }

    /* approximate marginal log Likelihood ................. */
    
    for (int i = 0; i < n; i++) {     
      
      switch (family) {
      case 1: /* cumulative link model */
	logLikCond_obs[i] = 
	  log(olmm_GLink(etaCLM[n + i] + etaRanefCLM[n + i], link) - 
	      olmm_GLink(etaCLM[i] + etaRanefCLM[i], link));
	logLikCond_sbj[subject[i]-1] += logLikCond_obs[i];
	break;
      case 2: case 3: /* baseline-category logit model Likelihood */
	for (int j = 0; j < nEta; j++)
	  sumBL[i] += exp(eta[n * j + i] + etaRanef[n * j + i]);
	logLikCond_obs[i] = -log(1.0 + sumBL[i]);
	if (y[i] < J)
	  logLikCond_obs[i] += eta[n * (y[i] - 1) + i] + 
	    etaRanef[n * (y[i] - 1) + i];
	logLikCond_sbj[subject[i]-1] += logLikCond_obs[i];
	break;
      }
    }
      
    for (int i = 0; i < N; i++) {
      logLik_sbj[i] += exp(logLikCond_sbj[i] + gq_weight);
    }
    
    if (numGrad == 0) {
	
      /* approximate score function ........................ */
      
      for (int i = 0; i < n; i++) {

	for (int j = 0; j < nEta; j++)
	  etaTmp[j] = eta[n * j + i] + etaRanef[n * j + i];

	/* hack to avoid numeric problems with DBL_MIN, DBL_MAX */
	logLikCond_modified = 
	  fabs(exp(logLikCond_obs[i])) < 1E-6 & gq_weight < 1E-6 ?
	  DBL_MAX : logLikCond_obs[i];
	
	/* calculate the Kronecker product of u_i and w_it */
	for (int j = 0; j < q; j++) {
	  for (int l = 0; l < nEta; l++) {
	    for (int m = 0; m < qEtaVar; m++) {
	      vecRanefTerm[q * j + l * qEtaVar + m] = 
		gq_nodes[j] * W[n * m + i];
	    }
	  }
	  for (int l = 0; l < qEtaInv; l++) {
	    vecRanefTerm[q * j + nEta * qEtaVar + l] = 
	      gq_nodes[j] * W[n * (qEtaVar + l) + i];
	  }
	}
	
	/* get vRanefTerm */
	F77_CALL(dgemv)("N", &lenVRanefCholFac, &lenVecRanefCholFac, 
			&one, ranefElMat, 
			&lenVRanefCholFac, vecRanefTerm, &i1, 
			&zero, vRanefTerm, &i1);
	
	/* fixed effect scores ........................ */
	
	/* predictor-variable fixed effects scores .... */
      
	switch (family) {
	case 1: 
	  if (y[i] < J) { /* score on upper concerned effects */
	    scoreCondVar[2] = 
	      olmm_gLink(etaCLM[n + i] + etaRanefCLM[n + i], link) / 
	      exp(logLikCond_modified);

	    for (int j = 0; j < pEtaVar; j++)
	      scoreCond_obs[n * (pEtaVar*(y[i]-1) + j) + i] +=
		scoreCondVar[2] * X[n * j + i];
	  }	    
	  if (y[i] > 1) { /* score on lower effects */
	    scoreCondVar[1] = 
	      -olmm_gLink(etaCLM[i] + etaRanefCLM[i], link) / 
	      exp(logLikCond_modified);

	    for (int j = 0; j < pEtaVar; j++) {
	      scoreCond_obs[n * (pEtaVar * (y[i]-2) + j) + i] +=
		scoreCondVar[1] * X[n * j + i];
	    }
	  }
	  break;
	case 2: case 3:
	  
	  for (int j = 0; j < nEta; j++) {
	    scoreCondVar[j] = 
	      (y[i]-1 == j ? 1.0 : 0.0) -
	      exp(etaTmp[j]) / (1.0 + sumBL[i]);

	    for (int l = 0; l < pEtaVar; l++) {
	      scoreCond_obs[n * (pEtaVar * j + l) + i] = 
		scoreCondVar[j] * X[n * l + i];
	    }
	  }
	  break;
	}

	/* predictor-invariant fixed effect scores .... */
	
	switch (family) {
	case 1:
	  scoreCondInv = /* will be used also for random effect scores */
	    (olmm_gLink(etaCLM[i + n] + etaRanefCLM[i + n], link) -
	     olmm_gLink(etaCLM[i] + etaRanefCLM[i], link)) / 
	    exp(logLikCond_modified);

	  for (int j = 0; j < pEtaInv; j++) {
	    scoreCond_obs[n * (nEta * pEtaVar + j) + i] =
	      scoreCondInv * X[n * (pEtaVar + j) + i];
	  }
	  break;
	case 2: /* baseline-category model */
	  scoreCondInv = 
	    (y[i] < J ? 1.0 : 0.0) -
	    sumBL[i] / (1.0 + sumBL[i]);

	  for (int j = 0; j < pEtaInv; j++) {
	    scoreCond_obs[n * (nEta * pEtaVar + j) + i] =
	      scoreCondInv * X[n * (pEtaVar + j) + i];	      
	  }
	  break;
	case 3: /* adjacent-categories model */
	  scoreCondInv = 
	    (y[i] < J ? 1.0 : 0.0) -
	      sumBL[i] / (1.0 + sumBL[i]);

	  for (int j = 0; j < pEtaInv; j++) {
	    for (int l = 0; l < nEta; l++) {
	      scoreCond_obs[n * (nEta * pEtaVar + j) + i] +=
		scoreCondVar[l] * (J - l - 1) * X[n * (j + pEtaVar) + i];
	    }
	  }
	  break;
	}
    
	/* random effect scores ....................... */

	for (int j = 0; j < q; j++) { /* columns of RanefCholFac */
	  
	  for(int l = j; l < q; l++) { /* rows of RanefCholFac */
	    
	    subsTmp = j * q - j * (j + 1) /2 + l;

	    if (l < qEtaVar * nEta) {
	      
	      /* predictor-variable random effect scores */
	      
	      tmpJ = (int)(floor(l / qEtaVar));
	      
	      switch (family) {
	      case 1: 
		if (y[i] < J & y[i] == tmpJ) {
		  scoreCond_obs[n * (p + subsTmp) + i] +=
		    scoreCondVar[2] * vRanefTerm[subsTmp];
		}
		if (y[i] > 1 & y[i] == tmpJ) {
		  scoreCond_obs[n * (p + subsTmp) + i] += 
		    scoreCondVar[1] * vRanefTerm[subsTmp];
		}
		break;
	      case 2: case 3: 
		scoreCond_obs[n * (p + subsTmp) + i] = 
		  scoreCondVar[tmpJ] * vRanefTerm[subsTmp];
		break;
	      }

	    } else {
	      
	      /* predictor-invariant random effect scores */
	      
	      switch (family) {
	      case 1:
		scoreCond_obs[n * (p + subsTmp) + i] =
		  scoreCondInv * vRanefTerm[subsTmp];
		break;
	      case 2:
		scoreCond_obs[n * (p + subsTmp) + i] =
		  scoreCondInv * vRanefTerm[subsTmp];
		break;
	      case 3:
		for (int m = 0; m < nEta; m++) {
		  scoreCond_obs[n * (p + subsTmp) + i] += 
		    scoreCondVar[m] * (J - m - 1) * vRanefTerm[subsTmp];
		}
		break;
	      }
	    }
	  }
	}
      }
      
      /* add terms to score_obs */      
      for (int i = 0; i < n; i++)
	for (int j = 0; j < nPar; j++)	
	  score_obs[n * j + i] += scoreCond_obs[n * j + i] * 
	    exp(logLikCond_sbj[subject[i]-1] + gq_weight);

    }
 
    gq_weight = 1; /* reset integration weight */
  }

  Free(etaRanef);
  Free(logLikCond_obs);
  Free(logLikCond_sbj);
  if (family == 1) Free(etaCLM);
  if (family == 1) Free(etaRanefCLM);
  if (family == 2 | family == 3) Free(sumBL);
  Free(scoreCond_obs);
  Free(scoreCond_sbj);

  /* add up score_obs to score_sbj (before computing info matrix!) */

  if (numGrad == 0) {

    for (int i = 0; i < n; i++)
      for (int j = 0; j < nPar; j++)
	score_sbj[N * j + subject[i]-1] += score_obs[n * j + i];
  }

  /* compute information matrix ................ */
  
  if (numHess == 0) {
 
    AllocVal(info, nPar * nPar, zero); /* reset info matrix*/
    double *hDerVec = Alloca(nPar, double), hInv2;
    
    for (int i = 0; i < N; i++) {    
      hInv2 = - weights_sbj[i] / (logLik_sbj[i] * logLik_sbj[i]);
      for (int j = 0; j < nPar; j++) hDerVec[j] = score_sbj[i + N * j];
      F77_CALL(dgemm)("N", "T", &nPar, &nPar, &i1, &hInv2, hDerVec,
		      &nPar, hDerVec, &nPar, &one, info, &nPar);
    }
  }

  /* finish score function calculations ........ */
  
  if (numGrad == 0) {

    for (int i = 0; i < n; i++)
      for (int j = 0; j < nPar; j++)
	score_obs[n * j + i] *= 1 / logLik_sbj[subject[i]-1];

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < nPar; j++) {
	score_sbj[N * j + i] *= 1 / logLik_sbj[i];
	score[j] += weights_sbj[i] * score_sbj[N * j + i];
      }
    }
  }

  /* finish Likelihood function calculations ..... */

  for (int i = 0; i < N; i++) {
    logLik_sbj[i] = weights_sbj[i] * log(logLik_sbj[i]);
    logLik[0] += logLik_sbj[i];
  }

  UNPROTECT(2);
  return R_NilValue;
}

/**
 * ---------------------------------------------------------
 * Expected random effects given fixed effects.
 *    
 * @param  x a olmm object
 *
 * @return R_NilValue
 * ---------------------------------------------------------
 */

SEXP olmm_update_u(SEXP x) {

  /* get grouping factor slot */  
  SEXP subject_0 = GET_SLOT(x, olmm_subjectSym), subject_1;
  PROTECT(subject_1 = coerceVector(subject_0, INTSXP));

  /* get targed variable slot */
  SEXP y_0 = GET_SLOT(x, olmm_ySym), y_1;
  PROTECT(y_1 = coerceVector(y_0, INTSXP));
  R_CheckStack();

  /* integer valued slots and pointer to factor valued slots */
  int *y = INTEGER(y_1), *subject = INTEGER(subject_1),
    *dims = DIMS_SLOT(x);
  R_CheckStack();

  /* numeric valued objects */
  double *X = X_SLOT(x), *W = W_SLOT(x), *u = U_SLOT(x),
    *offset = OFFSET_SLOT(x), *eta = ETA_SLOT(x), 
    *fixef = FIXEF_SLOT(x), *ranefCholFac = RANEFCHOLFAC_SLOT(x), 
    *logLik_sbj = LOGLIKSBJ_SLOT(x), *logLik = LOGLIK_SLOT(x), 
    *info = INFO_SLOT(x),
    *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x);
  R_CheckStack();

  /* set constants (dimensions of vectors etc.) */
  const int n = dims[n_POS], N = dims[N_POS], 
    p = dims[p_POS], 
    pEtaVar = dims[pEtaVar_POS], pEtaInv = dims[pEtaInv_POS], 
    q = dims[q_POS], qEta = dims[qEta_POS], 
    qEtaVar = dims[qEtaVar_POS], qEtaInv = dims[qEtaInv_POS],
    J = dims[J_POS], nPar = dims[nPar_POS], nEta = dims[nEta_POS],
    nGHQ = dims[nGHQ_POS], nQP = dims[nQP_POS], 
    family = dims[fTyp_POS], link = dims[lTyp_POS], 
    verb = dims[verb_POS],
    lenVecRanefCholFac = q * q, lenVRanefCholFac = q * (q + 1) / 2;

  /* variables for matrix operations */
  int i1 = 1; 
  double one = 1.0, zero = 0.0, tmp;

  /* define internal vectors */
  double *etaCLM = (double*) NULL,
    *etaRanefCLM = (double*) NULL,
    *sumBL = (double*) NULL,
    *etaRanef = Calloc(n * nEta, double),
    *ranefVec = Alloca(q, double);
  R_CheckStack();

  AllocVal(etaRanef, n * nEta, zero);

  switch (family) {
  case 1: 
    etaCLM = Calloc(n * 2 , double);
    etaRanefCLM = Calloc(n * 2 , double);
    for (int i = 0; i < n; i++) {
      etaCLM[i] /* lower */
	= y[i] > 1 ? eta[n * (y[i] - 2) + i] : -DBL_MAX; 
      etaCLM[i + n] = /* upper */
	y[i] < J ? eta[n * (y[i] - 1) + i] : DBL_MAX;
    }
    break;
  case 2: case 3:
    sumBL = Calloc(n, double);
    break;
  }
  R_CheckStack();

  AllocVal(u, N * q, zero);

  /* Gauss-Hermite quadrature */

  double *gq_nodes = Alloca(q, double), gq_weight = 1,
    *ranef = Alloca(q, double),
    *logLikCond_obs = Calloc(n, double),
    *logLikCond_sbj = Calloc(N, double);

  for (int k = 0; k < nQP; k++) {
    
    /* clear temporary objects */
    AllocVal(logLikCond_obs, n, zero);
    AllocVal(logLikCond_sbj, N, zero);
    if (family == 2 | family == 3) AllocVal(sumBL, n, zero);

    /* prepare integration nodes and weights */    
    for (int i = 0; i < q; i++) {
      gq_nodes[i] = ghx[nQP * i + k];
      gq_weight *= ghw[nQP * i + k];
    }
    gq_weight = log(gq_weight);

    /* multiply ranefCholFac with actual nodes */
    F77_CALL(dgemv)("N", &q, &q, &one, ranefCholFac, &q, 
		    gq_nodes, &i1, &zero, ranefVec, &i1);

    /* ranefVec (vector) to ranef (matrix) */
    for (int i = 0; i < nEta; i++) {
      Memcpy(ranef + i * qEta, ranefVec + i * qEtaVar, qEtaVar);
      Memcpy(ranef + i*qEta + qEtaVar, ranefVec + qEtaVar*nEta, qEtaInv);
    }
    
    /* set predictor-invariant random effects of adjacent-category model 
       to use the Likelihood of the baseline-category model */
    if (family == 3) {
      for (int i = 0; i < qEtaInv; i++) {
	for (int j = 0; j < nEta; j++) {
	  ranef[j * (qEtaVar + qEtaInv) + qEtaVar + i] *=
	    J - j - 1; 
	}
      }
    }

    /* compute contribution of random effects to linear predictor */
    F77_CALL(dgemm)("N", "N", &n, &nEta, &qEta, &one, W, &n, 
		    ranef, &qEta, &zero, etaRanef, &n);
 
    /* use etaRanefCLM to accelerate computations */
    if (family == 1) {
      for (int i = 0; i < n; i++) {
	etaRanefCLM[i] /* lower */
	  = y[i] > 1 ? etaRanef[n * (y[i] - 2) + i] : -DBL_MAX; 
	etaRanefCLM[i + n] = /* upper */
	  y[i] < J ? etaRanef[n * (y[i] - 1) + i] : DBL_MAX;
      }
    }
        
    for (int i = 0; i < n; i++) {     

      /* approximate log-Likelihood */

      switch (family) {
      case 1: /* cumulative link model */
	logLikCond_obs[i] = 
	  log(olmm_GLink(etaCLM[i + n] + etaRanefCLM[i + n], link) - 
	      olmm_GLink(etaCLM[i] + etaRanefCLM[i], link));
	logLikCond_sbj[subject[i]-1] += logLikCond_obs[i];
	break;
      case 2: case 3: /* baseline-category logit model Likelihood */
	for (int j = 0; j < nEta; j++)
	  sumBL[i] += exp(eta[i + j * n] + etaRanef[i + j * n]);
	logLikCond_obs[i] = -log(1.0 + sumBL[i]);
	if (y[i] < J)
	  logLikCond_obs[i] += eta[i + n * (y[i] - 1)] + 
	    etaRanef[i + n * (y[i] - 1)];
	logLikCond_sbj[subject[i]-1] += logLikCond_obs[i];
	break;
      }
    }
    
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < q; j++) {
	u[i + j * N] += gq_nodes[j] * exp(logLikCond_sbj[i] + gq_weight);
      }
    }

    gq_weight = 1; /* reset integration weight */
  }

  /* divide u by marginal Likelihood ............. */

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < q; j++) {
      u[i + j * N] /= exp(logLik_sbj[i]);
    }
  }

  Free(etaRanef);
  Free(logLikCond_obs);
  Free(logLikCond_sbj);
  if (family == 1) Free(etaCLM);
  if (family == 1) Free(etaRanefCLM);
  if (family == 2 | family == 3) Free(sumBL);

  UNPROTECT(2);
  return R_NilValue;
}


/**
 * ---------------------------------------------------------
 * Marginal predicition.
 *    
 * @param  x a olmm object
 * @param  p a matrix

 * @return R_NilValue
 * ---------------------------------------------------------
 */

SEXP olmm_pred_marg(SEXP x, SEXP eta, SEXP W, SEXP n, SEXP pred) {

  double *reta = REAL(eta), *rW = REAL(W), *rpred = REAL(pred);

  const int rn = INTEGER(n)[0];

  /* integer valued slots and pointer to factor valued slots */
  int *dims = DIMS_SLOT(x);
  R_CheckStack();

  /* numeric valued objects */
  double *ranefCholFac = RANEFCHOLFAC_SLOT(x),
    *ghw = GHW_SLOT(x), *ghx = GHX_SLOT(x);
  R_CheckStack();

  /* set constants (dimensions of vectors etc.) */
  const int q = dims[q_POS], qEta = dims[qEta_POS], 
    qEtaVar = dims[qEtaVar_POS], qEtaInv = dims[qEtaInv_POS],
    J = dims[J_POS], 
    nEta = dims[nEta_POS], nGHQ = dims[nGHQ_POS], nQP = dims[nQP_POS], 
    family = dims[fTyp_POS], link = dims[lTyp_POS], 
    verb = dims[verb_POS];

  /* variables for matrix operations etc. */
  int i1 = 1;
  double one = 1.0, zero = 0.0, 
    gq_weight = 1, sumBL = 0.0;

  /* define internal objects */
  double *etaRanef = Calloc(rn * nEta, double),
    *gq_nodes = Alloca(q, double),
    *predCond = Calloc(rn * J, double),
    *ranefVec = Alloca(q, double),
    *ranef = Alloca(qEta * nEta, double);
  R_CheckStack();
  
  AllocVal(rpred, rn * J, zero);

  for (int k = 0; k < nQP; k++) {

    AllocVal(predCond, rn * J, zero);

    /* prepare integration nodes and weights */    
    for (int i = 0; i < q; i++) {
      gq_nodes[i] = ghx[nQP * i + k];
      gq_weight *= ghw[nQP * i + k];
    }
    gq_weight = log(gq_weight);

    /* multiply ranefCholFac with actual nodes */
    F77_CALL(dgemv)("N", &q, &q, &one, ranefCholFac, &q, 
		    gq_nodes, &i1, &zero, ranefVec, &i1);

    /* ranefVec (vector) to ranef (matrix) */
    for (int i = 0; i < nEta; i++) {
      Memcpy(ranef + i * qEta, ranefVec + i * qEtaVar, qEtaVar);
      Memcpy(ranef + i*qEta + qEtaVar, ranefVec + qEtaVar*nEta, qEtaInv);
    }

    /* set predictor-invariant random effects of adjacent-category models 
       to use the Likelihood of the baseline-category model */
    if (family == 3) {
      for (int i = 0; i < qEtaInv; i++) {
	for (int j = 0; j < nEta; j++) {
	  ranef[(qEtaVar + qEtaInv) * j + qEtaVar + i] *=
	    J - j - 1; 
	}
      }
    }

    /* compute contribution of random effects to linear predictor */
    F77_CALL(dgemm)("N", "N", &rn, &nEta, &qEta, &one, rW, &rn, 
		    ranef, &qEta, &zero, etaRanef, &rn);

    /* approximate marginal probabilities .................. */
    
    for (int i = 0; i < rn; i++) { 
      for (int j = 0; j < J; j++) {

	switch (family) {
	case 1: /* cumulative link model */
	    predCond[rn * j + i] = 
	      log(olmm_GLink(j < (J - 1) ? reta[rn * j + i] + etaRanef[rn * j + i] : DBL_MAX, link) - olmm_GLink(j > 0 ? reta[rn * (j - 1) + i] + etaRanef[rn * (j - 1) + i] : -DBL_MAX, link));
	  break;
	case 2: case 3: /* baseline-category logit model Likelihood */
	  sumBL = 0;
	  for (int l = 0; l < nEta; l++)
	    sumBL += exp(reta[rn * l + i] + etaRanef[rn * l + i]);
	  predCond[rn * j + i] = -log(1.0 + sumBL);
	  if (j < (J - 1))
	    predCond[rn * j + i] += reta[rn * j + i] + etaRanef[rn * j + i];
	  break;
	}
      }
    }
    
    for (int i = 0; i < rn; i++) {
      for (int j = 0; j < J; j++) {
	rpred[rn * j + i] += exp(predCond[rn * j + i] + gq_weight);
      }
    }
    
    gq_weight = 1; /* reset integration weight */
  }
  Free(etaRanef);
  Free(predCond);
  return R_NilValue;
}
