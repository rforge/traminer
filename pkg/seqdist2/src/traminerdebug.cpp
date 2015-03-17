#include "TraMineR.h"
int TRAMINER_DEBUG_LEVEL=TRAMINER_DEBUG_LEVEL_DEFAULT;

extern "C"{
	SEXP setSDDebugLevel(SEXP level){
		TRAMINER_DEBUG_LEVEL= INTEGER(level)[0];
		return R_NilValue;
	}
	
	SEXP getSDDebugLevel(){
		return ScalarReal(TRAMINER_DEBUG_LEVEL);
	}
}

