#include "TWEDdistance.h"
TWEDdistance::TWEDdistance(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS)
	:OMdistance(normS, Ssequences, seqdim, lenS), 
	nu(0), lambda(0){
}
TWEDdistance::TWEDdistance(TWEDdistance *dc)
	:OMdistance(dc), 
	nu(dc->nu), lambda(dc->lambda){
}

void TWEDdistance::setParameters(SEXP params){
	OMdistance::setParameters(params);

	nu = REAL(getListElement(params, "nu"))[0];
	lambda = REAL(getListElement(params, "lambda"))[0];
}

TWEDdistance::~TWEDdistance(){
}
double TWEDdistance::distance(const int&is, const int& js){
double minimum=0, i_warp=0, j_warp=0, sub=0;//, lenmax=0;
    //etats compar�s
    int i_state, j_state;
    int i_m1_state, j_m1_state; // BH Need previous state also
    double cost, maxpossiblecost;
    int i=1;
    int j=1;
    int m=slen[is]+1;
    int n=slen[js]+1;
    int prefix=0;

    // BH Jun  2 2013 23:56:32: Stripping common prefixes is not appropriate for TWED
    // This is probably because it looks one token back, so perhaps stripping up 
    // to prefix-1 would work, TODO
    //    printf("Dealing with common prefix\n");

    // while (i<m&&j<n&&sequences[MINDICE(is,i-1,nseq)]==sequences[MINDICE(js,j-1,nseq)]) {
    //     i++;
    //     j++;
    //     prefix++;
    // }
    //+1 pour correspondre ? la matrice F

    //    printf("Dealing with FMAT\n");

    while (i<m) {
        j=prefix+1;
        //        printf("i loop: %d %d\n",i,j);
        while (j<n) {
          //          printf("j loop: %d %d\n",i,j);
            i_state=sequences[MINDICE(is,i-1,nseq)];
            j_state=sequences[MINDICE(js,j-1,nseq)];
            // printf("States: %d %d\n",i_state,j_state);
            if (i==1) { // BH: set previous to dummy value if before start, else previous
              i_m1_state=1;
                } else {
              i_m1_state=sequences[MINDICE(is,i-2,nseq)];
            }
            // printf("States i/j and - 1: %d %d %d %d\n",i_state,j_state,i_m1_state,j_m1_state);
            if (j==1) {
              j_m1_state=1;
                } else {
              j_m1_state=sequences[MINDICE(js,j-2,nseq)];
            }
            // printf("test 1\n");
//            j_m1_state=sequences[MINDICE(js,j-2,nseq)]; // needs a test for i==1 | j==1
            if ((i_state == j_state) && (i_m1_state == j_m1_state)) {
                cost = 0;
            } else {
                cost = scost[MINDICE(i_state,j_state,alphasize)]
                     + scost[MINDICE(i_m1_state,j_m1_state,alphasize)];
                //      				Rprintf("costs = %d %d, %d => %f \n",MINDICE(i_state,j_state,alphasize),i_state,j_state,cost);
            }

            i_warp = fmat[MINDICE(i-prefix,j-1-prefix,fmatsize)] + 
              scost[MINDICE(j_state,j_m1_state,alphasize)]
              + nu + lambda;

            j_warp = fmat[MINDICE(i-1-prefix,j-prefix,fmatsize)]+
              scost[MINDICE(i_state,i_m1_state,alphasize)]
              + nu + lambda;

            sub = fmat[MINDICE(i-1-prefix,j-1-prefix,fmatsize)]+ cost + 2*nu*abs(i-j);
            //            printf("i j iw ij match: %d %d %5.2f %5.2f %5.2f\n",i,j,i_warp, j_warp, sub);

            minimum = i_warp;
            if (j_warp < minimum) minimum = j_warp;

            if (sub < minimum) minimum = sub;

            fmat[MINDICE(i-prefix,j-prefix,fmatsize)]=minimum;
            j++;
        }
        i++;
    }//Fmat build

    // for (i=1; i<m; i++) {
    //   printf("%c", 65+sequences[MINDICE(is,i-1,nseq)]);
    // }
    // printf("\n");
    // for (i=1; i<m; i++) {
    //   printf("%c", 65+sequences[MINDICE(js,i-1,nseq)]);
    // }
    // printf("\n");
    // for (i=0; i<m; i++) {
    //   for (j=0; j<n; j++) {
    //     printf(" %5.2f", fmat[MINDICE(i,j,fmatsize)]);
    //   }
    //   printf("\n");
    // }


    //    printf("Finished with FMAT\n");

    m--;
    n--;
    //Warning! m and n decreased!!!!!
    maxpossiblecost=abs(n-m)*lambda+maxscost*fmin2((double)m,(double)n); // BH: indel replaced with lambda, probably incorrect May 30 2013 15:57:59
    return normalizeDistance(fmat[MINDICE(m-prefix,n-prefix,fmatsize)], maxpossiblecost, m, n);
}

