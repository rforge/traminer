#include "eventseq.h"
#include <sstream>

using namespace std;

/** Sequence finalizer, used by R to free memory
*/
void finalizeSequence(SEXP ptr) {

    Sequence *s;
    ASSIGN_TMRSEQ_TYPE(s,ptr);
    delete s;
}
/** CLASS SequenceEventNode
  Contain one event in an indiviudal sequence, reponsible for deleting next event in sequence
*/

//recursive add
void SequenceEventNode::addEvent(const int &e,const double &t) {
    if (this->hasNext()){
    	if(this->next->greaterThan(e, t-this->gap)){
			SequenceEventNode * s= new SequenceEventNode(e,t-this->gap);
			this->next->gap-=t-this->gap;
			s->setNext(this->next);
			this->next=s;
    	}else this->next->addEvent(e,t-this->gap);
    }
    else {
        this->next=new SequenceEventNode(e,t-this->gap);
    }
}

void Sequence::addEvent(const int &e,const double &t) {
    if (this->hasEvent()){
    	if(this->event->greaterThan(e,t)){
			this->event->setGap(this->event->getGap()-t);
    		SequenceEventNode * s=new SequenceEventNode(e,t);
    		s->setNext(this->event);
    		this->event=s;
    	}else{
			this->event->addEvent(e,t);
    	}
    } else {
        this->event=new SequenceEventNode(e,t);
    }
}
SequenceEventNode * SequenceEventNode::copy() {
    SequenceEventNode *s=new SequenceEventNode(this->type, this->gap);
    if (this->hasNext())s->next=this->next->copy();
    return s;
}
Sequence * Sequence::copy() {
    Sequence *s=new Sequence(this->idpers,this->dict);
    if (this->hasEvent())s->event=this->event->copy();
    return s;
}
///CLASS Sequence
///Represent an individual sequence
///CTor
Sequence::Sequence(const int&id, EventDictionary* ed):dict(ed), obsTime(-1), weight(1) { //personnal time 0, type =0 (root)
//    this->ns=NULL;
	this->dict->addSequence();
    this->idpers=id;
    this->event=NULL;
}
///Dtor clear all events and next sequences.
Sequence::~Sequence() {
    //if (ns!=NULL)delete ns;
    if (this->event!=NULL) delete event;
    this->dict->removeSequence();
    if(this->dict->shouldDelete())delete this->dict;
}

string Sequence::sprint() {
	ostringstream oss;
	oss.precision(2);
    //if(!this->isGeneric())n = sprintf(buffer, (char*)"[%i] ",this->idpers);
    //Rprintf((char*)"Current buffer %s\n",buffer);
    if (this->hasEvent()) {
        this->event->sprint(oss, true, !this->isGeneric(), (*this->dict), this->obsTime);
    }
    return oss.str();

}
void Sequence::print() {
    string r=this->sprint();
    //Rprintf((char *)"%s %i",buffer,r);
    REprintf((char *)"%s\n",r.c_str());
}
void SequenceEventNode::sprint(ostringstream &oss, const bool& start, const bool &printGap, const EventDictionary& ed, const double & remainingTime) {
    if (start) {
        if (this->gap>0&&printGap) {
			oss << this->gap << "-(" << ed.find(this->type)->second;
            //tmp=sprintf(&buffer[index],(char*)"%.2f-",this->gap);
        } else {
            //tmp=ed.sprint(&buffer[index],"(",this->type);
			oss << "(" << ed.find(this->type)->second;
        }


    } else if (this->gap>0) {
        if (printGap) {
			oss << ")-" << this->gap << "-(" << ed.find(this->type)->second;
        } else {
			oss << ")-(" << ed.find(this->type)->second;
        }
    } else {
		oss << "," << ed.find(this->type)->second;
    }
    if (this->hasNext()) {
        this->next->sprint(oss, false, printGap, ed, remainingTime-this->gap);
    } else {
    	if(remainingTime>0){
			oss << ")-" << (remainingTime-this->gap);
    	}
    	else{
    		oss << ")";
    	}

    }
}



double Sequence::first_occurence(Sequence * s, const double &maxGap, const double& windowSize, const double & ageMin, const double & ageMax, const double & ageMaxEnd) {
    if (!this->hasEvent()||!s->hasEvent()) return -1;
    double age=0;
    SequenceEventNode *sen=s->getEvent();
    while (sen!=NULL) {
        age+=sen->getGap();
        if (age>ageMax)return -1;
        if (age>=ageMin&&this->event->checkType(sen)&&this->event->count(sen,maxGap, windowSize,ageMaxEnd,0,age)>0) {
            return age;
        }
        sen=sen->getNext();
    }
    return -1;
}
int Sequence::count(Sequence * s, const double &maxGap, const double& windowSize, const double & ageMin, const double & ageMax, const double & ageMaxEnd) {
    if (!this->hasEvent()||!s->hasEvent()) return 0;
    SequenceEventNode * sen=s->getEvent();
    double age=0;
    int c=0;
    while (sen!=NULL) {
        age+=sen->getGap();
        if (age>ageMax)break;
        if (age>=ageMin&&this->event->checkType(sen)) {
            c+=this->event->count(sen, maxGap, windowSize,ageMaxEnd, 0,age);
        }
        sen=sen->getNext();
    }
    return c;
}
//Internal method, check, how many starting from this point already checked (we are looking for the rest of subsequence
int SequenceEventNode::count(SequenceEventNode * s, const double &maxGap, const double& windowSize, const double & ageMaxEnd, const double& gapConsumed, const double& currentAge) {
    int c=0;
    if (!this->hasNext())return 1;
    //Rprintf("Checked %i and %i %f\n", this->type,s->type, this->next->gap);
    SequenceEventNode * sen=s->getNext();
    if (this->next->gap==0) {
        //Rprintf("Equality %i and %i\n", this->next->type,sen->type);
        while (sen!=NULL&&sen->gap==0) {
            //Rprintf("Equality %i and %i\n", this->next->type,sen->type);
            if (this->next->checkType(sen)) {
                c+=this->next->count(sen, maxGap, windowSize, ageMaxEnd,gapConsumed,currentAge);
            }
            sen=sen->getNext();
        }
    } else {
        while (sen!=NULL&&sen->gap==0) {//Looking for next gap
            sen=sen->getNext();
            //Rprintf("Skipping %i\n", sen->type);
        }
        if (sen==NULL)return 0; //we didn't match
        double g=0;
        while (sen!=NULL) {
            //	Rprintf("Step %i and %i\n", this->next->type,sen->type);
            g+=sen->gap;
            //Time constraints, we don't need to look deeper
            if (g>maxGap||(g+gapConsumed)>windowSize)return c;
            if ((currentAge+g)>ageMaxEnd)return c;
            if (this->next->checkType(sen)) {
                c+=this->next->count(sen, maxGap, windowSize, ageMaxEnd,g+gapConsumed,g+currentAge);
            }
            sen=sen->getNext();
        }
    }
    return c;
}
bool Sequence::contain(const EventSet& es, const bool& exclude){
	if (!this->hasEvent()) return false;
    SequenceEventNode * sen=this->getEvent();
    while (sen!=NULL) {
		if(es.contain(sen->getType())){
			if(!exclude) return true;
		}else if(exclude)return false;
        sen=sen->getNext();
    }
    return exclude;
}
/*bool Sequence::notContain(const EventSet& es){
	if (!this->hasEvent()) return true;
    SequenceEventNode * sen=this->getEvent();
    while (sen!=NULL) {
		if(es.contain(sen->getType()))return false;
        sen=sen->getNext();
    }
    return true;
}*/
