#include "NMSMSTSoftdistanceII.h"
NMSMSTSoftdistanceII::NMSMSTSoftdistanceII(SEXP normS, SEXP Ssequences, SEXP seqdim, SEXP lenS):SUBSEQdistance(normS, Ssequences, seqdim, lenS){
	rowsize=maxlen+1;
	int matrixsize=rowsize*rowsize;
	e1=new double[matrixsize];
	e=new double[matrixsize];
}
NMSMSTSoftdistanceII::NMSMSTSoftdistanceII(NMSMSTSoftdistanceII *dc):SUBSEQdistance(dc), softmatch(dc->softmatch), alphasize(dc->alphasize){
	rowsize=maxlen+1;
	int matrixsize=rowsize*rowsize;
	e1=new double[matrixsize];
	e=new double[matrixsize];
}
NMSMSTSoftdistanceII::~NMSMSTSoftdistanceII(){
	delete [] e1;
	delete [] e;

}
void NMSMSTSoftdistanceII::setParameters(SEXP params){
	this->softmatch=REAL(getListElement(params, "softmatch"));
	this->alphasize=INTEGER(getListElement(params, "alphasize"))[0];
	SUBSEQdistance::setParameters(params);
}
void NMSMSTSoftdistanceII::computeattr(const int&is, const int& js){
	int mrows=slen[is]+1; //length of seq is
	int ncols=slen[js]+1;//length of seq js
	int m=mrows-1,n=ncols-1, ij=0;
	double s=0.0, t=0.0, phi_k=0.0;
	int seqi=-1, seqj=-1; //state in seq is pos i and in seq js, pos j
	double sf = 0.0;
	int k=0;
	
	/// Initialization of the algorithm
	for (int i=0;i<m;i++) {
		seqi=sequences[MINDICE(is,i,nseq)]; //state in seq is, pos i
		for (int j=0;j<n;j++) {
			ij = MINDICE(i,j,rowsize); // position ij in a matrix of size rowsize
			seqj = sequences[MINDICE(js,j,nseq)]; // state in seq js, pos j
			//if(seqi==sequences[MINDICE(js,j,nseq)]) {
			sf =softmatch[MINDICE(seqi,seqj,alphasize)]; //softmatching coefficient between seqj and seqi
			e1[ij]=sf; //initialize the matrix using the softmatching coefficient
			e[ij]=sf;
			s+=sf; //total of minimum shared time
			if(s==DBL_MAX){//ensure we do not have numerical problems
				error(" [!] Number of subsequences is getting too big"); 
			} 
		}
	}
	/// Initialization of the algorithm
	/// Last column and last row of the matrix =0.  
	for(int i=0; i<m;i++){
		ij = MINDICE(i,ncols-1,rowsize);
		e1[ij]=0.0;
		e[ij]=0.0;
	}
	for(int i=0; i<ncols;i++) {
		ij = MINDICE(mrows-1, i,rowsize);
		e1[ij]=0.0;
		e[ij]=0.0;
	}
	TMRLOG(5,"\n\ne1\n");
	TMRLOGMATRIX(5, e1, mrows,ncols, rowsize);
	TMRLOG(5,"\n\ne\n");
	TMRLOGMATRIX(5, e, mrows,ncols, rowsize);
	
	///End of initialization of the algorithm
	// Store the value for k=0
	this->kvect[k]=s+1;
	if(s==0.0){
		return;
	}
	while(mrows !=0 && ncols != 0) {
		k++;
		// sum up lower part of matrix e (matching) and t (minimum shared time) to update the matrix
		for (int irow=0;irow<mrows;irow++) {
			s = 0.0;
			for (int jcol=ncols-1;jcol>=0;jcol--) {
				ij = MINDICE(irow,jcol,rowsize);
				t = s; 
				s += e[ij];
				e[ij] = t;
			}
		}
		phi_k =0;
		for (int jcol=0;jcol<ncols;jcol++) {
			s = 0;
			for (int irow=mrows-1;irow>=0;irow--) {
				ij = MINDICE(irow, jcol, rowsize);
				t = s; 
				s += e[ij];
				e[ij] = t * e1[ij];
				phi_k += e[ij];
			}
		}
		if(phi_k==0){
			return;
		}
		this->kvect[k] = phi_k;
		if(phi_k == DBL_MAX){ // Ensure we do not have numerical errors
			error(" [!] Number of subsequences is getting too big");
		} 
		mrows--;
		ncols--;
		TMRLOG(5,"\n\ne\n");
		TMRLOGMATRIX(5, e, mrows,ncols, rowsize);
	}
}


	/***********************************************
	* Original version of the algorithm from Chesa 
	* Cees Elzinga
	* Calculates the minimal shared time (MST) spent
	* in all matching subsequences from strings 
	* seq[iseq] and seq[jseq]
	***********************************************/
/* static R_INLINE double NMSMSTdistance(int* slen,const int &is,const int &js, const int&nseq, int* sequences, const int &maxlen, 
								double * hmat, double * vmat, double * zmat, double * selfmatvect, double *ijvect ) {

	int mrows=slen[is]+1;
	int ncols=slen[js]+1;
	int m=mrows-1,n=ncols-1;
	int[,] e1=new int[mrows,ncols];
	int[,] e=new int[mrows,ncols];
	double[,] t1=new double[mrows,ncols];
	double[,] t=new double[mrows,ncols]; 
	double st=0.0,sum=0.0,temp=0.0,temps=0.0;
	double sumt,tempt,temptt;
	int seqi=-1;
	double ti;
	for (int i=0;i<m;i++) {
		seqi=DataSpec.seq[iseq][i];
		ti=DataSpec.quant[iseq][i];
		for (int j=0;j<n;j++) {
			if(seqi==DataSpec.seq[jseq][j]) {
				e1[i,j]=1;
				e[i,j]=1;
				t1[i,j]=Math.Min(ti,DataSpec.quant[jseq][j]);
				t[i,j]=t1[i,j];
				st+=t[i,j];
			}
		}
	}
	if(st==0.0){
		return st;
	}
	while(mrows !=0 && ncols != 0) {
		for (int irow=0;irow<mrows;irow++) {
			sum=0.0;
			sumt=0.0;
			for (int jcol=ncols-1;jcol>=0;jcol--) {
				temp=sum;
				tempt=sumt;
				sum+=e[irow,jcol];
				sumt+=t[irow,jcol];
				e[irow,jcol]=(int)temp;
				t[irow,jcol]=tempt;
			}
		}
		temps=0.0;
		temptt=0.0;
		for (int jcol=0;jcol<ncols;jcol++) {
			sum=0;
			sumt=0.0;
			for (int irow=mrows-1;irow>=0;irow--) {
				temp=sum;
				tempt=sumt;
				sum+=e[irow,jcol];
				sumt+=t[irow,jcol];
				e[irow,jcol]=(int)temp*e1[irow,jcol];
				t[irow,jcol]=e1[irow,jcol]*(e[irow,jcol]*t1[irow,jcol]+tempt);
				temps+=e[irow,jcol];
				temptt+=t[irow,jcol];
			}
		}
		if(temps==0.0){
			return st;
		}
		st+=temptt;
		mrows--;
		ncols--;
	} 
	return st;
}
 */
 