#include<R.h>
#include "eventseq.h"
#include "prefixtree.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "eventdictionary.h"
#include "constraint.h"
#include <Rmath.h>
#include "TraMineR.h"
#include <string>

/**
	tmrsequence build a sequence obect and return an external pointer to that object
*/




extern "C" {

    /**Forward declarations of functions, See below for more explanations*/
    /*SEXP tmrsequence(SEXP idpers, SEXP time, SEXP event, SEXP classname);
    SEXP tmrsequenceseveral(SEXP idpers, SEXP time, SEXP event,SEXP endEvent,
							SEXP classname, SEXP dictionnary);
    SEXP tmrsequencestring(SEXP seq);
    SEXP tmrsequencegetid(SEXP seq);
    SEXP tmrfindsubsequences(SEXP seqs,SEXP maxGap, SEXP windowSize,
                             SEXP ageMinBegin, SEXP ageMaxBegin, SEXP ageMaxEnd,
                             SEXP minSupport, SEXP maxSubseqSize, SEXP classname);
    SEXP tmrmatrixsubseqinseq(SEXP subseqs, SEXP seqs,SEXP maxGap, SEXP windowSize,
                              SEXP ageMinBegin, SEXP ageMaxBegin,SEXP ageMaxEnd,
                              SEXP countMethod);
    //Exported functions
    static R_CallMethodDef TMRSEQUENCE_CallDefs[] = {
        {"tmrsequence", (DL_FUNC) tmrsequence, 4},
        {"tmrsequenceseveral", (DL_FUNC) tmrsequenceseveral, 5},
        //  {"tmrprintsequence", (DL_FUNC) tmrprintsequence, 1},
        {"tmrsequencestring", (DL_FUNC) tmrsequencestring, 1},
        {"tmrsequencegetid", (DL_FUNC) tmrsequencegetid, 1},
        {"tmrfindsubsequences", (DL_FUNC) tmrfindsubsequences, 9},
        {"tmrmatrixsubseqinseq", (DL_FUNC) tmrmatrixsubseqinseq, 8},
        {NULL}
    };*/
    void R_init_TraMineR(DllInfo *info) {
        //   TMRSEQUENCE_type_tag= install("TMRSEQUENCE_TYPE_TAG");
        //R_registerRoutines(info, NULL, TMRSEQUENCE_CallDefs, NULL, 0);
    }

    /**
    	Build one sequence obect, a given idpers, time should be double and event integer
    	Ascending order is expect for time and event (if time is equal)
    */
    SEXP tmrsequence(SEXP idpers, SEXP time, SEXP event, SEXP classname, SEXP seq) {
        Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
		EventDictionary * ed=s->getDictionary();
		//Get pointers
        double * t=REAL(time);
        int *ev=INTEGER(event);
        int len=length(time);
        if (len!=length(event))error("Time and event vector arent of the same size");
        int id=INTEGER(idpers)[0], i;
        if (len==0)return R_NilValue;
        //Build sequence
        s=new Sequence(id,ed);
        //For each pair (time,event)
        for (i=0;i<len;i++) {
            //add the event
            s->addEvent(ev[i],t[i]);
        }
        //Return the Sequence as R object
        return makeTMRSequence(s,classname);
    }

//        return ans;
    /**
         	Build several sequences obects, a given idpers, time should be double and event integer
         	Ascending order is expect for time and event (if time is equal) grouped by idpers
         */
    SEXP tmrsequenceseveral(SEXP idpers, SEXP time, SEXP event, SEXP endEvent,SEXP classname, SEXP dictionnary) {
    	//Create the dictionnary
    	EventDictionary * ed= new EventDictionary(dictionnary);
    	bool obsTime=!isNull(endEvent);
    	int eEvent=0;
    	if(obsTime){
    		eEvent=INTEGER(endEvent)[0];
    	}
        //Time pointer
        double * t=REAL(time);
        //events and ids
        int *ev=INTEGER(event),*ids=INTEGER(idpers);
        //lengthes
        int totlen=length(time);
        //should all be the same size
        if (totlen!=length(event)||totlen!=length(idpers))error("Time ,idpers and event vector should have the same size");
        if (totlen==0)return R_NilValue;
        int id=ids[0], i,idlen=1,idpos=0;
        int lastID=id;
        //Counting number of distinct ID's
        for (i=0;i<totlen;i++) {
            if (ids[i]!=lastID) {
                lastID=ids[i];
                idlen++;
            }
        }
        lastID=id;
        if (idlen<1)error("Not enough sequences");
        //Rprintf((char*)"totlen %i : idlen: %i\n", totlen,idlen);
        SEXP ans,tmpseq;
        //List to return
        PROTECT(ans=allocVector(VECSXP, idlen));
        //building first sequence
        Sequence *s=new Sequence(id, ed);
        //For each pair (time,event)
        for (i=0;i<totlen;i++) {
            id=ids[i]; //Get current ID
            if (id!=lastID) { //If new, store old sequence
                tmpseq = makeTMRSequence(s, classname);
                SET_VECTOR_ELT(ans,idpos,tmpseq);//Put in vector
                idpos++;//position in vector
                s=new Sequence(id,ed);//Build new sequence
                lastID=id;
            }
            if(obsTime&&ev[i]==eEvent){
				s->setObsTime(t[i]);
            } else {
				s->addEvent(ev[i],t[i]);
            }
        }
        //Store last built sequence
        tmpseq = makeTMRSequence(s,classname);
        SET_VECTOR_ELT(ans,idpos,tmpseq);
        //Unprotect vector ans
        UNPROTECT(1);
        return ans;//Return ans
    }
	SEXP tmrsequencecontainevent(SEXP seqs, SEXP eventList, SEXP exclude) {

		EventSet es;
		es.add(eventList);
        //events and ids
        int numseq=length(seqs);
        bool excl=INTEGER(exclude)[0]==1;
        SEXP seq, ret;
        PROTECT(ret = allocVector(LGLSXP, numseq));
        int *pret=LOGICAL(ret);
        Sequence *s =NULL;
        //lengthes
		for (int i=0;i<numseq;i++) {
                seq=VECTOR_ELT(seqs,i);
                ASSIGN_TMRSEQ_TYPE(s,seq);
                pret[i]=s->contain(es, excl);
                //Rprintf((char*)"Added %i seq, node=%i\n",i,TreeEventNode::getNodeCount());
		}
		UNPROTECT(1);
        return ret;
    }

    SEXP tmrsequencegetid(SEXP seq) {
        Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
        return ScalarInteger(s->getIDpers());
    }
    SEXP tmrsequencegetlength(SEXP seq) {
        Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
        return ScalarReal(s->getObsTime());
    }
	SEXP tmrsequencegetweight(SEXP seq) {
        Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
        return ScalarReal(s->getWeight());
    }
    SEXP tmrsequencesetlength(SEXP seqs, SEXP time) {
    	 double * t=REAL(time);
        //events and ids
        int numseq=length(seqs);
        SEXP seq;
        Sequence *s =NULL;
        //lengthes
        if(length(time)!=numseq)error("Time and seq vector should have the same size");
		for (int i=0;i<numseq;i++) {
                seq=VECTOR_ELT(seqs,i);
                ASSIGN_TMRSEQ_TYPE(s,seq);
                s->setObsTime(t[i]);
                //Rprintf((char*)"Added %i seq, node=%i\n",i,TreeEventNode::getNodeCount());
		}
        return R_NilValue;
    }
	SEXP tmrsequencesetweight(SEXP seqs, SEXP weight) {
    	 double * w=REAL(weight);
        //events and ids
        int numseq=length(seqs);
        SEXP seq;
        Sequence *s =NULL;
        //lengthes
        if (length(weight)!=numseq) error("Weight and seq vector should have the same size");
		for (int i=0;i<numseq;i++) {
                seq=VECTOR_ELT(seqs,i);
                ASSIGN_TMRSEQ_TYPE(s,seq);
                s->setWeight(w[i]);
                //Rprintf((char*)"Added %i seq, node=%i\n",i,TreeEventNode::getNodeCount());
		}
        return R_NilValue;
    }
    SEXP tmrsequencegetdictionary(SEXP seq){
		Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
        return s->getDictionary()->getDictionary();
    }
/**Return a string representation of a sequence*/
    SEXP tmrsequencestringinternal(SEXP seq) {
        Sequence *s =NULL;
        ASSIGN_TMRSEQ_TYPE(s,seq);
        std::string buffer = s->sprint();
        return mkChar(buffer.c_str());
    }
    /**Return a string representation of a sequence*/
    SEXP tmrsequencestring(SEXP seq) {
        SEXP str=tmrsequencestringinternal(seq);
        SEXP ret;
        PROTECT(ret = allocVector(STRSXP, 1));
        SET_STRING_ELT(ret, 0, str);
        UNPROTECT(1);
        return ret;
    }


    /**Main function find frequent subsequences*/
    SEXP tmrfindsubsequences(SEXP seqs,SEXP maxGap, SEXP windowSize,SEXP ageMinBegin, SEXP ageMaxBegin, SEXP ageMaxEnd, SEXP countMethod, SEXP minSupport, SEXP maxSubseqSize, SEXP classname) {
        //Initializing parameters

    	Constraint * cst = new Constraint(REAL(maxGap)[0],REAL(windowSize)[0],REAL(ageMinBegin)[0],REAL(ageMaxBegin)[0],REAL(ageMaxEnd)[0], REAL(countMethod)[0]);
    	//REprintf((char*)"branches/nico 1\n\n");

     // MOVED TO CONSTRAINT OBJ
     //  double wSize=REAL(windowSize)[0],mGap=REAL(maxGap)[0];
     //  double aMin=REAL(ageMinBegin)[0],aMax=REAL(ageMaxBegin)[0],aMaxEnd=REAL(ageMaxEnd)[0];
	/*	double wSize = cst->getwindowSize();
		double mGap = cst->getmaxGap();
		//double mGap=REAL(maxGap)[0];
		double aMin=cst->getageMinBegin();
		double aMax=cst->getageMaxBegin();
		double aMaxEnd=cst->getageMaxEnd();
    	REprintf((char*)"wSize = %f\n", wSize);
    	REprintf((char*)"maxgap = %f\n", mGap);
    	REprintf((char*)"aMin = %f\n", aMin);
    	REprintf((char*)"aMax = %f\n", aMax);
    	REprintf((char*)"aMaxEnd = %f\n", aMaxEnd);
	 */
        int maxK=INTEGER(maxSubseqSize)[0],k=1;
		double mSupport=REAL(minSupport)[0];

        //Default values implies no limit (actually biggest possible limit)
        // MOVED TO CONSTRAINT OBJ
       // if (wSize==-1)wSize=DBL_MAX;
      //  if (mGap==-1)mGap=DBL_MAX;
      //  if (aMax==-1)aMax=DBL_MAX;
       // if (aMaxEnd==-1)aMaxEnd=DBL_MAX;
        if (maxK==-1)maxK=INT_MAX;
        SEXP seq;
        int numseq=length(seqs);
        Sequence * s=NULL;
        int lastNodeCount;
        PrefixTree * root= new PrefixTree();
        EventDictionary * ed=NULL;
        //Adding one event to subseq at a time
        for (k=1;k<=maxK;k++) {
            //Clear support stored in tree
            root->clearSupport();
            lastNodeCount=TreeEventNode::getNodeCount();
            //REprintf((char*)"Step %i:\n     Adding sequences (size: %i)\n",k,TreeEventNode::getNodeCount());
            //Add every sequence to the tree
            for (int i=0;i<numseq;i++) {
                seq=VECTOR_ELT(seqs,i);
                ASSIGN_TMRSEQ_TYPE(s,seq);
                if(ed==NULL)ed=s->getDictionary();
                //root->addSequence(s,mGap,wSize,aMin,aMax,aMaxEnd,k);
                root->addSequence(s, cst, k);
              //  Rprintf((char*)"Added %i seq, node=%i\n",i,TreeEventNode::getNodeCount());
            }
            //REprintf((char*)"     Simplifying tree (size: %i)\n",TreeEventNode::getNodeCount());
            //root->print();
            //return ScalarLogical(TRUE);
            //Simplify tree
            root->simplifyTree(mSupport);
            //REprintf((char*)"     Tree simplified (size: %i [added: %i])\n",TreeEventNode::getNodeCount(), (TreeEventNode::getNodeCount()-lastNodeCount));
            if (TreeEventNode::getNodeCount()-lastNodeCount==0)break;
        }
        // root->print();
        //root->print();
        //Tree size (number of node=number of frequent subsequences)
        //REprintf((char*)"Counting subseq...");
        int returnsize=root->countSubsequence(mSupport);
        //Rprintf((char*)"Counting subseq (%i)\n",returnsize);
        SEXP ans, supp,subseq;
        PROTECT(ans=allocVector(VECSXP, 2)); //Allocate memory
        PROTECT(supp=allocVector(REALSXP, returnsize)); //Allocate memory
        PROTECT(subseq=allocVector(VECSXP, returnsize)); //Allocate memory
        int index=0;
        //REprintf((char*)"(%i)\nRetrieving subsequences...",returnsize);
        root->getSubsequences(subseq,REAL(supp),&index,classname, ed); //Extracting all subsequences
        //REprintf((char*)"OK\n");
        SET_VECTOR_ELT(ans,0,supp);
        SET_VECTOR_ELT(ans,1,subseq);
        UNPROTECT(3);
		delete cst;
        delete root;
        return ans;
    }



    /**Count the number of time we can find subseq in the sequence seq (with given maxGap and windowSize)
    	ageMinBegin, ageMaxBegin and ageMaxEnd permits to handle time constraints
    	countMethod define the kind of return:
    		1: number of occurences
    		2: presence-absence
    		3: age at first occurrence
    */
    SEXP tmrmatrixsubseqinseq(SEXP subseqs, SEXP seqs,SEXP maxGap, SEXP windowSize,SEXP ageMinBegin, SEXP ageMaxBegin,SEXP ageMaxEnd, SEXP countMethod) {
        double wSize=REAL(windowSize)[0],mGap=REAL(maxGap)[0];
        double aMin=REAL(ageMinBegin)[0],aMax=REAL(ageMaxBegin)[0], aMaxEnd=REAL(ageMaxEnd)[0];
        int cMethod=INTEGER(countMethod)[0];
        if (wSize==-1)wSize=DBL_MAX;
        if (mGap==-1)mGap=DBL_MAX;
        if (aMax==-1)aMax=DBL_MAX;
        if (aMaxEnd==-1)aMaxEnd=DBL_MAX;
        Sequence *s =NULL, *sub=NULL;
        int nsub=length(subseqs);
        int ns=length(seqs);
        SEXP ans;
        SEXP subseq,seq, namesubseq, nameseq,dimnames;
        PROTECT(ans = allocMatrix(REALSXP, ns, nsub));
        double *matrix=REAL(ans);
        PROTECT(namesubseq= allocVector(STRSXP, nsub));
        PROTECT(nameseq= allocVector(STRSXP, ns));
        for (int j=0;j<ns;j++) {
            seq=VECTOR_ELT(seqs,j);
            SET_STRING_ELT(nameseq, j,tmrsequencestringinternal(seq));
        }
        for (int i=0;i<nsub;i++) {
            subseq=VECTOR_ELT(subseqs,i);
            ASSIGN_TMRSEQ_TYPE(sub,subseq);
            SET_STRING_ELT(namesubseq, i,tmrsequencestringinternal(subseq));
            //Rprintf("Processing ");
            //sub->print();
            for (int j=0;j<ns;j++) {
                seq=VECTOR_ELT(seqs,j);
                ASSIGN_TMRSEQ_TYPE(s,seq);
                //Rprintf("Counting on ");
                //s->print();
                //int counting=sub->count(s,mGap,wSize,aMin,aMax);
                //Rprintf("Counted %i\n",counting);
                //matrix[j+i*ns]=counting;
                switch (cMethod) {
                case 1:
                    matrix[j+i*ns]=sub->count(s,mGap,wSize,aMin,aMax, aMaxEnd);
                    break;
                case 2:
                    matrix[j+i*ns]=fmin2(1.0,sub->count(s,mGap,wSize,aMin,aMax, aMaxEnd));
                    break;
                case 3:
                    matrix[j+i*ns]=sub->first_occurence(s,mGap,wSize, aMin, aMax, aMaxEnd);
                    break;
                default:
                    matrix[j+i*ns]=sub->count(s,mGap,wSize,aMin,aMax, aMaxEnd);
                    break;
                }

            }
        }
        PROTECT(dimnames = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0,nameseq);
        SET_VECTOR_ELT(dimnames, 1, namesubseq);
        setAttrib(ans, R_DimNamesSymbol, dimnames);
        UNPROTECT(4);
        return ans;
    }
	/**Find each times an events appears, return a matrix with ncol = maximum number of the specified event
    */
    SEXP tmreventinseq(SEXP seqs, SEXP Sevent) {
        int event=INTEGER(Sevent)[0];
        Sequence *s =NULL;
        int ns=length(seqs);
        SEXP ans;
        SEXP seq;
		int nseqevent=0, maxnevent=1;
		SequenceEventNode * sen=NULL;
		//Start by counting the maximum number of time, the event appears in a given sequences
		for (int i=0;i<ns;i++) {
            seq=VECTOR_ELT(seqs,i);
			ASSIGN_TMRSEQ_TYPE(s,seq);
            if(s->hasEvent()){
				sen=s->getEvent();
				nseqevent=0;
				while(sen!=NULL){
					if(sen->getType()==event){
						nseqevent++;
					}
					sen=sen->getNext();
				}
				if(nseqevent>maxnevent)maxnevent=nseqevent;
            }
        }
		TMRLOG(4, "Maximum numbers of event %d is %d", event, maxnevent);
        PROTECT(ans = allocMatrix(REALSXP, ns, maxnevent));
        double *matrix=REAL(ans);
		double age=0;
		//looking up for events ages
        for (int i=0;i<ns;i++) {
            seq=VECTOR_ELT(seqs,i);
			ASSIGN_TMRSEQ_TYPE(s,seq);
			nseqevent=0;
            if(s->hasEvent()){
				sen=s->getEvent();
				age=0;
              
				while(sen!=NULL){
					age += sen->getGap();
					if(sen->getType()==event){
						matrix[MINDICE(i,nseqevent,ns)]=age;
						nseqevent++;
					}
					sen=sen->getNext();
				}
            }
			//Filling non used values with -1
			while(nseqevent<maxnevent){
				matrix[MINDICE(i,nseqevent,ns)]=-1;
				nseqevent++;
			}
        }
        UNPROTECT(1);
        return ans;
    }

}