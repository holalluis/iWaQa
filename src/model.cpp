/*
 *  model.cpp
 *  generic model functionality
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 
 
#include "model.h"

#include "mathutils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <float.h>

//flexible data conversion through a stringstream
template<typename to, typename from> to lexical_cast(from const &x) 
{
  std::stringstream os;
  to ret;

  os << x;
  os >> ret;

  return ret;  
}

//-------------------------------------------------------------------------------------------------------------

//Utility to compare the end of strings (needed in parameter aliasing in case of @0

bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

//-------------------------------------------------------------------------------------------------------------

//Utility to replace parts of a string
void searchAndReplace(std::string& value, std::string const& search,std::string const& replace)
{
    std::string::size_type  next;

    for(next = value.find(search);        // Try and find the first match
        next != std::string::npos;        // next is npos if nothing was found
        next = value.find(search,next)    // search for the next match starting after
                                          // the last match that was found.
       )
    {
        // Inside the loop. So we found a match.
        value.replace(next,search.length(),replace);   // Do the replacement.
        next += replace.length();                      // Move to just after the replace
                                                       // This is the point were we start
                                                       // the next search from. 
    }
}

//-------------------------------------------------------------------------------------------------------------

#pragma mark Parameter manager

iWQParameterManager::iWQParameterManager()
{
	//
}

//-------------------------------------------------------------------------------------------------------------

iWQParameterManager::~iWQParameterManager()
{
	//notify (and kick out) bound clients
	for(int i=0; i<mBoundClients.size(); i++){
		iWQModel * act_client=mBoundClients[i];
		if(act_client){
			act_client->detach();
		}
	}
	//delete storage	
	for(int i=0; i<mLocalParams.size(); i++){
		delete mLocalParams[i];
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::defineParam(std::string key)
{
	//allocate among mLocalParams 
	double * d=new double;
	mLocalParams.push_back(d);
	*d=0.0;
	//load into the name storage
	mParams[key]=d;
	//init the distribution strorage
	mOrderedDistributions.push_back(NULL);	//no distribution by default
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::setValueForParam(double value, std::string key)
{
	double * par=mParams[key];
	if(par){
		*par=value;
		for(int i=0; i<mBoundClients.size(); i++){
			mBoundClients[i]->updateParameters();
		}
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::initParam(std::string key, double value)
{
	defineParam(key);
	setValueForParam(value,key);
}

//-------------------------------------------------------------------------------------------------------------

double iWQParameterManager::valueForParam(std::string key)
{
	double * par=mParams[key];
	if(par){
		return *par;
	}
	return 0.0;
}

//-------------------------------------------------------------------------------------------------------------

bool iWQParameterManager::hasValueForParam(std::string key)
{
	return (mParams.find(key)!=mParams.end());
}

//-------------------------------------------------------------------------------------------------------------

int iWQParameterManager::numberOfParams()
{
	//only local parameters
	return mLocalParams.size();
}

//-------------------------------------------------------------------------------------------------------------

std::vector<double> iWQParameterManager::plainValues()
{
	std::vector<double> result;
	for(int i=0; i<mLocalParams.size(); i++){
		result.push_back(*(mLocalParams[i]));
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

std::vector<std::string> iWQParameterManager::namesForPlainValues()
{
	std::vector<std::string> result;
	for(int i=0; i<mLocalParams.size(); i++){
		double * ptr=mLocalParams[i];
		std::map<std::string, double *>::iterator it;
		std::string act_name="< UNDEFINED >";
		for(it=mParams.begin(); it!=mParams.end(); ++it){
			if(it->second==ptr){
				act_name=makeLiteralFlaggedStr(it->first);
				break;
			}
		}
		result.push_back(act_name);
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::setPlainValues(std::vector<double> new_values)
{
	for(int i=0; i<mLocalParams.size() && i<new_values.size(); i++){
		*(mLocalParams[i])=new_values[i];
	}
	for(int i=0; i<mBoundClients.size(); i++){
		mBoundClients[i]->updateParameters();
	}
}	

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::setPlainValues(double * new_values, int numvals)
{
	for(int i=0; i<mLocalParams.size() && i<numvals; i++){
		*(mLocalParams[i])=new_values[i];
	}
	for(int i=0; i<mBoundClients.size(); i++){
		mBoundClients[i]->updateParameters();
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::bindRequest(iWQModel * client)
{
	//insert client to client-list
	//1. check if it is there
	for(int i=0; i<mBoundClients.size(); i++){
		if(mBoundClients[i]==client){
			return;
		}
	}
	mBoundClients.push_back(client);
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::detachRequest(iWQModel * client)
{
	//remove client from client-list
	for(int i=0; i<mBoundClients.size(); i++){
		if(mBoundClients[i]==client){
			mBoundClients.erase(mBoundClients.begin()+i);
			return;
		}
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::initFromFile(std::string filename)
{
	//read keys and values from the file
	FILE * in;
	char Buffer [512];
		
	in=fopen(filename.c_str(),"rt");
	if(in==NULL){
		return;
	}

	//clean previous values: not necessarily
	//mLocalParams.clear();
	//mParams.clear();
	
	std::map<std::string, double> new_values;	//temporary container for new values
	
	//get new values from the file
	while(fgets(Buffer,512,in)>0){
	
		std::string line=Buffer;
		//filter out whitespace
		std::remove(line.begin(), line.end(), ' ');
		std::remove(line.begin(), line.end(), '\t');
		
		//check comment position
		size_t commentpos=line.find('#');
		
		//truncate to non-comment
		if(commentpos!=std::string::npos){
			line=line.substr(0,commentpos);
		}
		
		//divide at the ':' 
		if(line.length()>0){
			size_t dividerpos=line.find(':');
			if(dividerpos!=std::string::npos){
				std::string parname=line.substr(0,dividerpos);
				std::string svalue=line.substr(dividerpos+1);
				const char * p1;
				char * p2;
				p1=svalue.c_str(); 
				double value=strtod(p1, &p2);
				if(p1!=p2){
					std::string keystr=makeFlaggedStr(parname);
					
					new_values[keystr]=value;
					//initParam(keystr,value);
				}
			}
		}
	}
	fclose(in);
	
	//now apply the changes
	std::map<std::string, double>::iterator it;
	for(it=new_values.begin(); it!=new_values.end(); ++it){
		//decide if this is a new parameter or not
		if(hasValueForParam(it->first)){
			setValueForParam(it->second, it->first);	//old: only update to keep the pointer bindings alive				
		}
		else{
			initParam(it->first, it->second);			//new: define from scratch
		}
		//remove those which did not appear in new_values: bad idea
		//because it would break existing pointer links and cause an access violation
		//when old links are resolved
	}
	
	//force refresh
	for(int i=0; i<mBoundClients.size(); i++){
		mBoundClients[i]->updateParameters();
	}
	
	return;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::initFromTabDelimitedFile(std::string filename)
{
	//read keys and values from the file
	FILE * in;
	char Buffer [512];
		
	in=fopen(filename.c_str(),"rt");
	if(in==NULL){
		return;
	}

	//clean previous values
	mLocalParams.clear();
	mParams.clear();
	
	//get new values from the file
	while(fgets(Buffer,512,in)>0){
	
		std::string line=Buffer;
		
		//check comment position
		size_t commentpos=line.find('#');
		
		//truncate to non-comment
		if(commentpos!=std::string::npos){
			line=line.substr(0,commentpos);
		}
		
		//scan the remaining
		double value;
		if(sscanf(line.c_str(),"%s\t%lf",Buffer,&value)==2){
			std::string parname=Buffer;
			std::string keystr=makeFlaggedStr(parname);
			initParam(keystr,value);
		}
	}
	fclose(in);
	
	//force refresh
	for(int i=0; i<mBoundClients.size(); i++){
		mBoundClients[i]->updateParameters();
	}
	
	return;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::saveToFile(std::string filename, bool tabdelimited)
{
	FILE * ofile=fopen(filename.c_str(),"wt");
	if(ofile==NULL){
		return;
	}
	std::map<std::string, double *>::iterator it;
	const char * pattern=(tabdelimited?"%s\t%g\n":"%s: %g\n");
	for(it=mParams.begin(); it!=mParams.end(); ++it){
		std::string act_key=it->first;
		double * act_dest=it->second;
		//try to disassemble the name
		if(act_key.length()>0 && act_key[0]=='#'){
			act_key = makeLiteralFlaggedStr(act_key);
		}
		if(act_dest){
			fprintf(ofile,pattern,act_key.c_str(),*act_dest);
		}
		else{
			printf("[Error]: Parameter %s has no value ptr.\n", act_key.c_str());
		}
	}
	fclose(ofile);
	return;
}

//-------------------------------------------------------------------------------------------------------------

std::string iWQParameterManager::makeLiteralFlaggedStr(std::string flaggedstr)
{
	std::string rootkey="";
	std::string flag="";
	decomposeFlaggedStr(flaggedstr, rootkey, flag);
	if(rootkey.length() && flag.length()){
		return rootkey + "[" + flag + "]";
	}
	else{
		return flaggedstr;
	}
}

//-------------------------------------------------------------------------------------------------------------

std::string iWQParameterManager::makeFlaggedStr(std::string literalflaggedstr)
{
	//filter out whitespace
	std::remove(literalflaggedstr.begin(), literalflaggedstr.end(), ' ');
	std::remove(literalflaggedstr.begin(), literalflaggedstr.end(), '\t');
	size_t closepos=literalflaggedstr.rfind(']');
	size_t openpos=literalflaggedstr.rfind('[');
	if(closepos!=std::string::npos && openpos!=std::string::npos && openpos>0 && openpos<closepos && closepos==literalflaggedstr.length()-1){
		std::string rootkey=literalflaggedstr.substr(0,openpos);
		std::string flag=literalflaggedstr.substr(openpos+1, closepos-openpos-1);
		if(rootkey.length()==0){
			return literalflaggedstr;
		}
		if(flag.length()==0){
			return rootkey;
		}
		else{
			return makeFlaggedStr(rootkey, flag);
		}
	}
	return literalflaggedstr;
}

//-------------------------------------------------------------------------------------------------------------
	
void iWQParameterManager::setValueForParam(double value, std::string key, std::string flag)
{
	std::string flagged=makeFlaggedStr(key, flag);
	setValueForParam(value, flagged);
}

//-------------------------------------------------------------------------------------------------------------

double iWQParameterManager::valueForParam(std::string key, std::string flag)
{
	std::string flagged=makeFlaggedStr(key, flag);
	return valueForParam(flagged);
}

//-------------------------------------------------------------------------------------------------------------

bool iWQParameterManager::hasValueForParam(std::string key, std::string flag)
{
	std::string flagged=makeFlaggedStr(key, flag);
	return hasValueForParam(flagged);
}

//-------------------------------------------------------------------------------------------------------------

std::string iWQParameterManager::makeFlaggedStr(std::string key, std::string flag)
{
	if(flag.size()==0){
		return key;
	}
	
	std::ostringstream os;
	os<<"# ";
	os<<key.length();
	os<<" ";
	os<<flag.length();
	os<<" @";
	std::string descriptorstr=os.str();
	
	return descriptorstr + key + flag;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::decomposeFlaggedStr(std::string flaggedstr, std::string & key, std::string & flag)
{
	key="";
	flag="";
	char firstchar='\0';
	int keylen=0;
	int flaglen=0;
	std::istringstream iss (flaggedstr);
	iss >> firstchar;
	//is it a valid composite?
	if(firstchar!='#'){
		key=flaggedstr;
		return;
	}
	iss >> keylen;
	iss >> flaglen;
	size_t atpos=flaggedstr.find('@');
	if(atpos==std::string::npos){
		//could not find the end of descriptor
		key=flaggedstr;
		return;
	}
	std::string content=flaggedstr.substr(atpos+1);
	if(content.length()!=keylen+flaglen){
		//sizes do not match
		key=flaggedstr;
		return;
	}
	key=content.substr(0,keylen);
	flag=content.substr(keylen,flaglen);
	
	return;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::defineParam(std::string key, std::string flag)
{
	std::string flagged=makeFlaggedStr(key, flag);
	defineParam(flagged);
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::initParam(std::string key, std::string flag, double value)
{
	std::string flagged=makeFlaggedStr(key, flag); 
	initParam(flagged, value);
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::setLimitsForParam(iWQLimits lim, std::string key, std::string domain)
{
	std::string flagged=makeFlaggedStr(key, domain);
	mLimits[flagged]=lim;
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::clearLimitsForParam(std::string key)
{
	mLimits.erase(key);
}

//-------------------------------------------------------------------------------------------------------------

bool iWQParameterManager::hasLimitsForParam(std::string key)
{
	
	return (mLimits.find(key)!=mLimits.end());
}

//-------------------------------------------------------------------------------------------------------------

iWQLimits iWQParameterManager::limitsForParam(std::string key)
{
	return mLimits[key];
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::linkDistributionToParam(iWQRandomGenerator * dist, std::string key, std::string flag)
{
	std::string fullkey=makeFlaggedStr(key, flag);
	if(dist && mParams.find(fullkey)!=mParams.end()){	//a valid distribution for a valid key
		mLinkedDistributions[fullkey]=dist;
		
		//load also to the ordered set
		int index=indexOfParam(fullkey);
		if(index>=0 && index<mOrderedDistributions.size()){
			mOrderedDistributions[index]=dist;
		}		
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::detachDistributionFromParam(std::string key, std::string flag)
{
	std::string fullkey=makeFlaggedStr(key, flag);
	mLinkedDistributions.erase(fullkey);
	int index=indexOfParam(fullkey);
	if(index>=0 && index<mOrderedDistributions.size()){
		mOrderedDistributions[index]=NULL;
	}
}

//-------------------------------------------------------------------------------------------------------------

iWQRandomGenerator * iWQParameterManager::distributionForParam(std::string key, std::string flag)
{
	std::string fullkey=makeFlaggedStr(key, flag);
	std::map<std::string, iWQRandomGenerator *>::iterator it=mLinkedDistributions.find(fullkey);
	if(it!=mLinkedDistributions.end()){
		return it->second;
	}
	else{
		return NULL;
	}
}

//-------------------------------------------------------------------------------------------------------------

double iWQParameterManager::logLikelihood(bool report)	
{
	//take all log likelihoods from the local parameters
	//std::map<std::string, double *>::iterator it;
	double loglikeli=0.0;
	//Optimized version
	int n=mLocalParams.size();
	for(int i=0; i<n; i++){
		iWQRandomGenerator * act_dist=mOrderedDistributions[i];
		if(act_dist){
			double act_likeli=act_dist->logLikeli(*mLocalParams[i]);
			if(report && (act_likeli== -DBL_MAX || isinf(act_likeli) || isnan(act_likeli))){
				//something is not likely at all...
				double * param_ptr = mLocalParams[i];
				std::map<std::string, double *>::iterator it;
				std::string parname="<Unknown parameter>";
				for(it=mParams.begin(); it!=mParams.end(); ++it){
					if(it->second==param_ptr){
						parname=it->first;
					}
				}
				printf("Prior log likelihood for %s=%lf is NA\n", parname.c_str(), *mLocalParams[i]);
			}
			
			//filter for infinitely unlikeliness
			if(act_likeli == -DBL_MAX){
				return -DBL_MAX;	//no reason to see the others
			}
			else{
				loglikeli+=act_likeli;
			}
		}
	}
	return loglikeli;
}

//-------------------------------------------------------------------------------------------------------------

double iWQParameterManager::logLikelihoodOfSet(std::vector<double> values)
{
	//need to revise this too (TODO)
	double loglikeli=0.0;
	std::vector<std::string> names=namesForPlainValues();
	for(int i=0; i<names.size() && i<values.size(); i++){
		std::string act_key=names[i];
		//ask if we have a generator for that name
		std::map<std::string, iWQRandomGenerator *>::iterator it2=mLinkedDistributions.find(act_key);
		if(it2!=mLinkedDistributions.end()){
			iWQRandomGenerator * act_dist=it2->second;
			if(act_dist){
				double act_likeli=act_dist->logLikeli(values[i]);
				loglikeli+=act_likeli;
			}
		}
	}
	return loglikeli;
}

//-------------------------------------------------------------------------------------------------------------

int iWQParameterManager::indexOfParam(std::string key)	//internal storage index, not necessarily the same as in plainValues, key is already flagged if applicable
{
	std::map<std::string, double *>::iterator it=mParams.find(key);
	if(it!=mParams.end()){
		double * ptr=it->second;
		//now look for that pointer in the storage 
		for(int i=0; i<mLocalParams.size(); i++){
			if(ptr==mLocalParams[i]){
				return i;
			}
		}
	}
	return -1;
}

//-------------------------------------------------------------------------------------------------------------

std::map<int, double> iWQParameterManager::valuesForParam(std::string key, std::string flag)
{
	std::map<int, double> result;
	
	//look through all parameter names
	std::map<std::string, double*>::iterator it;
	
	std::string rootkey;
	std::string flag2;
		
	for(it = mParams.begin(); it != mParams.end(); ++it){
		std::string rawkey = it->first;
		decomposeFlaggedStr(rawkey, rootkey, flag2);
		//check if rootkey matches
		if(key.compare(rootkey)==0 && flag2.size()>=3){
			//split flag
			size_t atpos=flag2.find('=');
			if(atpos!=std::string::npos){
				std::string rootflag = flag2.substr(0,atpos);
				//check if rootflag matches
				if(rootflag.compare(flag)==0){
					std::string flagflag = flag2.substr(atpos+1);
					//cast it to int
					int i = lexical_cast<int, std::string>(flagflag);
					result[i] = *(it->second);
				}
			}
			
		}
	}
	
	return result;
	
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::clearAllParams()
{
	for(int x=0; x<mLocalParams.size(); x++){
		delete mLocalParams[x];
	}
	mLocalParams.clear();
	mParams.clear();
	mLinkedDistributions.clear();
	mLimits.clear();
	mOrderedDistributions.clear();
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::clearParam(std::string key)
{
	int index = indexOfParam(key);
	if(index>=0 && index<mLocalParams.size()){
		double * d=mLocalParams[index];
		mLocalParams.erase(mLocalParams.begin()+index);
		mParams.erase(key);
		mLinkedDistributions.erase(key);
		mLimits.erase(key);
		mOrderedDistributions.erase(mOrderedDistributions.begin()+index);
		delete d;
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQParameterManager::clearParam(std::string key, std::string flag)
{
	std::string flagged=makeFlaggedStr(key, flag);
	clearParam(flagged);
}

//#############################################################################################################

#pragma mark Diagnostic and working functions for iWQModel variables

double * iWQDelta(iWQModel * aModel, double * variable)
{
	return aModel->mDerivLocations.operator[](variable);
}

//-------------------------------------------------------------------------------------------------------------

double * iWQDiagnosticDelta(iWQModel * aModel, double * variable)
{
	if(aModel->mDerivLocations.find(variable)!=aModel->mDerivLocations.end()){
		//normal operation
		return aModel->mDerivLocations.operator[](variable);
	}
	else{
		//error happened
		printf("[Error]: Unknown variable referenced with d() or F() in %s.\n",aModel->mTypeId.c_str());
		aModel->mVarDeltaErrorIndicator=true;
		return aModel->mFooVarDeltaContainer;
	}
}

//-------------------------------------------------------------------------------------------------------------

#pragma mark Model

iWQModel::iWQModel(std::string type) : mTypeId (type) , mModelId ("<unnamed>")
{
	mParentParameterManager=NULL;
	
	//Init the built-in Runge Kutta Fehlberg solver:
	//Solver tables
	mA=new double [6];
	mB=new double * [6];
	for(int i=0; i<6; i++){
		mB[i]=new double [5];
	}
	mC=new double [6];
	mD=new double [6];
	mF=new double * [6];
	
	//Constant tables
	mA[1]=0.25;				mA[2]=0.375;			 mA[3]=12.0/13.0;			mA[4]=1.0;				mA[5]=0.5;
	mB[1][0]=0.25;
	mB[2][0]=3.0/32.0;		mB[2][1]=9.0/32.0;
	mB[3][0]=1932.0/2197.0;	mB[3][1]=-7200.0/2197.0; mB[3][2]=7296.0/2197.0;
	mB[4][0]=439.0/216.0;	mB[4][1]=-8.0;			 mB[4][2]=3680.0/513.0;		mB[4][3]=-845.0/4104.0;
	mB[5][0]=-8.0/27.0;		mB[5][1]=2.0;			 mB[5][2]=-3544.0/2565.0;	mB[5][3]=1859.0/4104.0;	mB[5][4]=-11.0/40.0;
	mC[0]=16.0/135.0;		mC[1]=0.0;				 mC[2]=6656.0/12825.0;		mC[3]=28561.0/56430.0;	mC[4]=-9.0/50.0;		mC[5]=2.0/55.0;
	mD[0]=1.0/360.0;		mD[1]=0.0;				 mD[2]=-128.0/4275.0;		mD[3]=-2197.0/75240.0;	mD[4]=1.0/50.0;			mD[5]=2.0/55.0;	
	
	//default execution type (not diagnostic)
	dFunction = iWQDelta;
	mFooVarDeltaContainer = new double;
}

//-------------------------------------------------------------------------------------------------------------

iWQModel::~iWQModel()
{
	std::map<double*,double*>::iterator it;
	for(it=mDerivLocations.begin(); it!=mDerivLocations.end(); ++it){
		delete it->second;
	}
	delete mFooVarDeltaContainer;
	
	//delete the RKF solver arrays
	delete [] mA;
	delete [] mC;
	delete [] mD;
	delete [] mF;
	for(int i=0; i<6; i++){
		delete [] mB[i];
	}
	delete [] mB;
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::defineParam(double * par, const char * parname)
{
	if(par){			
		//this is where model parameters are defined
		std::string parnamestr=parname;
		std::string toreplace = "__AT0";
		if(hasEnding(parnamestr,toreplace)){			//replace __AT0 ending with @0 in string name
			parnamestr = parnamestr.substr(0, parnamestr.size()-toreplace.size()) + "@0";
		}
		mParams[parnamestr]=par;
		mParamInitState[parnamestr]=false;
		*par=0.0;
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::setValueForParam(double value, std::string key)
{
	double * par=mParams[key];
	mParamInitState[key]=true;
	if(par){
		*par=value;
	}
}

//-------------------------------------------------------------------------------------------------------------

double iWQModel::valueForParam(std::string key)
{
	//no need to ask the parent, everything should be synchronised with that
	double * par=mParams[key];
	if(par){
		return *par;
	}
	
	return 0.0;
}

//-------------------------------------------------------------------------------------------------------------

bool iWQModel::hasValueForParam(std::string key)
{
	bool result=false;
	if(mParentParameterManager){
		result=mParentParameterManager->hasValueForParam(key);
	}
	if(!result){
		result=(mParams.find(key) != mParams.end());
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::bind(iWQParameterManager * par)
{
	if(par){
		mParentParameterManager=par;
		par->bindRequest(this);
		updateParameters();		//initiated locally to get the latest values
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::detach()
{
	if(mParentParameterManager){
		mParentParameterManager->detachRequest(this);
	}
	mParentParameterManager=NULL;
}

//-------------------------------------------------------------------------------------------------------------

bool iWQModel::isBound() const
{
	return (mParentParameterManager!=NULL);
}

//-------------------------------------------------------------------------------------------------------------

iWQParameterManager * iWQModel::sharedManager()
{
	return mParentParameterManager;
}	

//-------------------------------------------------------------------------------------------------------------

void iWQModel::updateParameters()
{
	//updates everything which is stored in mParameterManager, preferring the flagged variants 
	if(!mParentParameterManager){
		return;
	}
	std::map<std::string, double *>::iterator it;
	for(it=mParams.begin(); it!=mParams.end(); ++it){
		std::string act_key=it->first;
		double * act_dest=it->second;
		//try out flags
		bool foundflagged=false;
		for(int i=0; i<mModelFlags.size(); i++){
			std::string act_flag=mModelFlags[i];
			if(mParentParameterManager->hasValueForParam(act_key, act_flag)){
				double new_val=mParentParameterManager->valueForParam(act_key, act_flag);
				*act_dest=new_val;
				foundflagged=true;
				break;
			}	
		}
		if(!foundflagged){
			//ask for the general one
			if(mParentParameterManager->hasValueForParam(act_key)){
				double new_val=mParentParameterManager->valueForParam(act_key);
				*act_dest=new_val;
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::defineVariable(double * var, const char * varname, bool delta)
{
	mVarLocations.push_back(var);
	mVarNames.push_back(varname);
	*var=0.0;
	double * deriv=new double;
	*deriv=0.0;
	mDerivLocations[var]=deriv;
	mShouldTakeDelta.push_back(delta);
}

//-------------------------------------------------------------------------------------------------------------

//for direct descendants (they use references to variables)
double * iWQModel::d(double & var)
{
	return dFunction(this,&var);
}

//-------------------------------------------------------------------------------------------------------------

//for indirect descendants (they can only use pointers to variables)
double * iWQModel::D(double * var)
{
	return dFunction(this,var);
}

//-------------------------------------------------------------------------------------------------------------

//alias for d to use with boundary fluxes
double * iWQModel::F(double & var)
{
	return d(var);
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::readVariables(double * from, int length)
{
	//copy variable values to the convenience storage
	if(!from){
		return;
	}
	for(int i=0; i<length && i<mVarLocations.size(); i++){
		double * var=mVarLocations[i];
		*var=from[i];
	}
	for(int i=mVarLocations.size(); i<length; i++){
		*mVarLocations[i]=0.0;
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::copyDerivatives(double * dest, int length)
{
	//place derivatives from the strage into dest
	if(!dest){
		return;
	}
	for(int i=0; i<length && i<mVarLocations.size(); i++){
		double * var=mVarLocations[i];
		double * deriv=mDerivLocations[var];
		dest[i]=*deriv;
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::setInitialValues(iWQInitialValues * initvals)
{
    if(!initvals){
        return;
    }
    for(int i=0; i<mVarNames.size(); i++){
        std::string act_name=mVarNames[i];
        
		//does initvals have a value for that key?
        double val=0.0;
		if(initvals->hasValueForVariable(act_name, mModelId)){
			val=initvals->valueForVariable(act_name, mModelId);	
		}
		else if(initvals->hasDefaultValueForVariable(act_name)){
			val=initvals->defaultValueForVariable(act_name);
		}
		
        *(mVarLocations[i])=val;
    }
}

//-------------------------------------------------------------------------------------------------------------

void iWQModel::defineInput(double * inp, const char * inputname)
{
	mInputLocations.push_back(inp);
	*inp=0.0;
	mInputNames.push_back(inputname);
}

//-------------------------------------------------------------------------------------------------------------

iWQStrings iWQModel::outputDataHeaders() const
{
	iWQStrings result;
	for(int i=0; i<mVarNames.size(); i++){
		if(!mShouldTakeDelta[i]){
			result.push_back(mVarNames[i]);
		}
	}
	for(int i=0; i<mVarNames.size(); i++){
		if(mShouldTakeDelta[i]){
			result.push_back(mVarNames[i]);
		}
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

iWQStrings iWQModel::inputDataHeaders() const
{
	iWQStrings result;
	for(int i=0; i<mInputNames.size(); i++){
		result.push_back(mInputNames[i]);
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

iWQStrings iWQModel::parameters() const
{
	iWQStrings result;
	std::map<std::string, double *>::const_iterator it;
	for(it=mParams.begin(); it!=mParams.end(); ++it){
		if(it->second){
			result.push_back(it->first);
		}
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

int iWQModel::numVariables() const
{
	int result=0;
	for(int i=0; i<mShouldTakeDelta.size(); i++){
		if(!mShouldTakeDelta[i]){
			result++;
		}
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

double * iWQModel::rwoutlet(std::string name) const
{
	//look up input
	double * result=NULL;
	for(int x=0; x<mInputNames.size(); x++){
		if(mInputNames[x]==name){
			return mInputLocations[x];
		}
	}
	return result;
}

//-------------------------------------------------------------------------------------------------------------

const double * iWQModel::routlet(std::string name) const
{
	//look up variable / input / local parameter
	double * result=NULL;
	for(int x=0; x<mVarNames.size(); x++){
		if(mVarNames[x]==name){
			return mVarLocations[x];
		}
	}
	for(int x=0; x<mInputNames.size(); x++){
		if(mInputNames[x]==name){
			return mInputLocations[x];
		}
	}
	std::map<std::string, double *>::const_iterator it = mParams.find(name);
	if(it!=mParams.end()){
		result=it->second;
	}
	return result;
}

//---------------------------------------------------------------------------------------------------------------

iWQStrings iWQModel::modelFlags() const
{
	return mModelFlags;
}

//---------------------------------------------------------------------------------------------------------------

void iWQModel::setModelFlags(iWQStrings flags)
{
	mModelFlags=flags;
	//refresh parameter values if needed
	updateParameters(); 
}		

//---------------------------------------------------------------------------------------------------------------

void iWQModel::setModelFlag(std::string flag)
{
	mModelFlags.push_back(flag);
}

//---------------------------------------------------------------------------------------------------------------

void iWQModel::setModelId(std::string newid)
{
	mModelId=newid;
}

//---------------------------------------------------------------------------------------------------------------

std::string iWQModel::modelId() const
{
	return mModelId;
}

//---------------------------------------------------------------------------------------------------------------

std::string iWQModel::modelType() const
{
	return mTypeId;
}

//---------------------------------------------------------------------------------------------------------------

bool iWQModel::verify()
{
	//main diagnostic function
	mVarDeltaErrorIndicator=false;
		
	bool result=true;
	
	//diagnose variables
	dFunction=iWQDiagnosticDelta;
	std::map<double *, double *>::iterator it;
	for(it=mDerivLocations.begin(); it!=mDerivLocations.end(); ++it){
		*(it->second)=0.0;
	}
	modelFunction(0.0);
	//check if we got something invalid in the delta container
	bool problematic=false;
	for(it=mDerivLocations.begin(); it!=mDerivLocations.end(); ++it){
		double value=*(it->second);
		if(isnan(value) || isinf(value)){
			//something strange happened
			problematic=true;
			std::string varname="<unknown variable>";
			std::string errtype=(isnan(value)?"NaN":"infinity");
			double * var_address=it->first;
			for(int j=0; j<mVarLocations.size(); j++){
				if(mVarLocations[j]==var_address){
					varname=mVarNames[j];
				}
			}
			printf("[Warning]: %s produced invalid derivative or flux value (%s) for %s.\n",mTypeId.c_str(),errtype.c_str(),varname.c_str());
		}
	}
	if(problematic){
		//print out parameter values for the invalid results
		printf("\tParameter values for this case:\n");
		std::map<std::string, double *>::iterator it3;
		for(it3=mParams.begin(); it3!=mParams.end(); ++it3){
			std::string act_key=it3->first;
			double * act_dest=it3->second;
			printf("\t\t%s:\t",act_key.c_str());
			if(act_dest){
				printf("%lf",*act_dest);
			}
			else{
				printf("<undefined>");
			}
			printf("\n");
		}
	}
	
	//check if the diagnostic delta function set the error flag
	result=(mVarDeltaErrorIndicator!=false)?false:result;
	
	//check parameter initialization
	std::map<std::string, double *>::iterator it2;
	for(it2=mParams.begin(); it2!=mParams.end(); ++it2){
		bool gotvalue=false;
		std::string act_key=it2->first;
		double * act_dest=it2->second;
		//try out flags
		bool foundflagged=false;
		for(int i=0; i<mModelFlags.size(); i++){
			std::string act_flag=mModelFlags[i];
			if(mParentParameterManager->hasValueForParam(act_key, act_flag)){
				gotvalue=true;
				foundflagged=true;
				break;
			}	
		}
		if(!foundflagged){
			//ask for the general one
			if(mParentParameterManager->hasValueForParam(act_key)){
				gotvalue=true;
			}
		}
		if(!gotvalue){
			//check if that has a local value
			if(mParamInitState[act_key]){
				gotvalue=true;
			}
		}
		if(!gotvalue){
			printf("[Warning]: Parameter %s of %s is not initialized.\n",act_key.c_str(),mModelId.c_str());
		}
	}
	
	//check if model run perturbs input storages
	int numinputs=mInputLocations.size();
	std::vector<double> input_backup;
	input_backup.assign(numinputs,0);
	std::vector<double> inputTestValues;
	inputTestValues.assign(numinputs,0);
	std::vector<bool> display_flag;
	display_flag.assign(numinputs,false);
	iWQRandomNormalGenerator r;
	r.setStdev(10);
	//make a backup copy of parameter values
	std::map<std::string, double> param_backup;
	std::map<std::string, double *>::iterator iter;
	for(iter=mParams.begin(); iter!=mParams.end(); ++iter){
		std::string act_key=iter->first;
		double * ptr=iter->second;
		if(ptr){
			param_backup[act_key] = *ptr;
		}
	}
	//make a backup from input values
	for(int i=0; i<numinputs; i++){
		input_backup[i]=*(mInputLocations[i]);
	}
	//perform test
	for(int j=0; j<10; j++){
		//make random inputs
		for(int i=0; i<numinputs; i++){
			inputTestValues[i]=r.generate();
		}
		//load values
		for(int i=0; i<numinputs; i++){
			*(mInputLocations[i])=inputTestValues[i];
		}
		//create random parameters
		for(iter=mParams.begin(); iter!=mParams.end(); ++iter){
			std::string act_key=iter->first;
			double * ptr=iter->second;
			if(ptr){
				*ptr = r.generate();
			}
		}
		//run the model function
		modelFunction(0.0);
		//examine input storages
		for(int i=0; i<numinputs; i++){
			if(*(mInputLocations[i])!=inputTestValues[i]){
				if(!display_flag[i]){
					printf("[Error]: %s changes the value in input container %s.\n",mTypeId.c_str(),mInputNames[i].c_str());
					display_flag[i]=true;
				}
				result=false;
			}
		}
	}
	//restore original params
	for(iter=mParams.begin(); iter!=mParams.end(); ++iter){
		std::string act_key=iter->first;
		double * ptr=iter->second;
		if(ptr){
			*ptr = param_backup[act_key];
		}
		else{
			printf("[Error]: Could not restore %s in %s.\n",act_key.c_str(),mTypeId.c_str());
		}
	}
	//restore original inputs
	for(int i=0; i<numinputs; i++){
		*(mInputLocations[i]) = input_backup[i];
	}
	
	// END OF TESTS
	// restore normal operation
	dFunction=iWQDelta;
	
	return result;
}

//---------------------------------------------------------------------------------------------------------------

bool iWQModel::solve1Step(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps)
{
	//verify parameters
	if(!verifyParameters()){	
		return false;
	}	
	
	if(!isStatic()){
		//Dynamic models: decide which solver to use
		//return solve1StepRungeKuttaFehlberg(xvon, xbis, yvon, hmin, eps);
		return solve1StepLSODA(xvon, xbis, yvon, hmin, eps);
	}
	else{
		//Simple solution for static models (no integration, no initial values)
		unsigned int numVariables=mVarLocations.size();
		double * ys;
		ys=new double [numVariables];
		modelFunction(xbis);						//FUNCTION CALLED
		copyDerivatives(ys,numVariables);	//get back the derivatives
		//reload new values into variables
		for(int k=0; k<numVariables; k++){
			*(mVarLocations[k])=ys[k];
		}
		delete [] ys;
		return true;
	}
}

//---------------------------------------------------------------------------------------------------------------

bool iWQModel::solve1StepRungeKuttaFehlberg(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps)
{
	double x, xs, h, hmax;
	double * y;
	double * ys;
	double * yhut;
	//double * inputs;						//PRELIMINARY TRIAL: is it useless?
	double Err, MaxErr, GrossErr, HNeu;
	int i, j, k;
	bool validityflag=true;
	
	unsigned int numVariables=mVarLocations.size();
	//unsigned int numInputs=mInputLocations.size();
	
	if(xbis<=xvon){
		return true;
	}
	
	//Storages for variables
	y=new double [numVariables];
	ys=new double [numVariables];
	yhut=new double [numVariables];
	
	for(i=0; i<=5; i++){
		mF[i]=new double [numVariables];
	}  
	
	//Storage for inputs
	
	hmax=xbis-xvon;
	h=hmax;
	
	xs=xvon;
	
	//initial variable values: will not modify anything if yvon is NULL
	setInitialValues(yvon);
	
	//get values back from the model into the plain storage
	for(k=0; k<numVariables; k++){
		ys[k]=(!mShouldTakeDelta[k])?*(mVarLocations[k]):0.0;
	}
	
	do{		
		readVariables(ys,numVariables);			//send the var data from ys to the derived class
		modelFunction(xs);						//FUNCTION CALLED
		copyDerivatives(mF[0],numVariables);	//get back the derivatives
		
		do{
			for(i=0; i<=5; i++){
				x = xs + mA[i] * h;
				for(k=0; k<numVariables; k++){
					y[k]=0.0;
					for(j=0; j<i; j++){
						y[k] += mB[i][j] * mF[j][k];
					}
					y[k] = h * y[k] + ys[k];
				}
		
				readVariables(y,numVariables);
				modelFunction(x);				//FUNCTION CALLED
				copyDerivatives(mF[i],numVariables);
			}
			GrossErr=0.0;
			for(k=0; k<numVariables; k++){
				yhut[k]=0.0;
				Err=0.0;
				for(i=0; i<=5; i++){
					yhut[k]+=mC[i]*mF[i][k];
					Err+=mD[i]*mF[i][k];
				}
				yhut[k]=h*yhut[k]+ys[k];
				Err=h*fabs(Err);
				if(Err>GrossErr){
					GrossErr=Err;
				}
			}
			MaxErr=h*eps;
			HNeu=(GrossErr!=0.0)?0.9*h*pow(MaxErr/GrossErr,0.25):hmax;
			if(HNeu<hmin){
				HNeu=hmin;
				validityflag=false;
			}
			//HNeu=(HNeu<hmin)?hmin:HNeu;
			if(GrossErr>MaxErr){
				h=HNeu;
			}
		}while(GrossErr>MaxErr && h>hmin);
		
		//Values at the end of the interval
		for(k=0; k<numVariables; k++){
			ys[k]=yhut[k];
		}
		xs+=h;
		
		//Reached the end of the step: export results
		if(xs==xbis){		
			//reload new values into variables
			for(k=0; k<numVariables; k++){
				*(mVarLocations[k])=(mShouldTakeDelta[k])?ys[k]/(xbis-xvon):ys[k];
			}
			break;
		}
		        
		h = HNeu;
        
		//Adjust step length to the remaining interval
		if(xs + h > xbis){
			h = xbis - xs;
		}
	}while(xs<xbis);
	
	//cleanup
	delete [] y;
	delete [] ys;
	delete [] yhut;
	
	for(i=0; i<=5; i++){
		delete [] mF[i];
	}
	
	return validityflag;
}

//---------------------------------------------------------------------------------------------------------------

bool iWQModel::solve1StepLSODA(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps)
{
	bool validityflag=true;
	
	unsigned int numVariables=mVarLocations.size();
	unsigned int numInputs=mInputLocations.size();
	
	if(xbis<=xvon){
		return true;
	}
	
	//initial variable values: will not modify anything if yvon is NULL
	setInitialValues(yvon);
	
	// now should create and call the LSODA integrator
	iWQLSODAIntegrator * integrator=new iWQLSODAIntegrator();
	validityflag=integrator->solve1Step(this, xvon, xbis, eps);
	
	delete integrator;
		
	return validityflag;
}

//--------------------------------------------------------------------------------------------------

//BEGIN DANGEROUS LOW-LEVEL FUNCTIONS
void iWQModel::resetState()
{
	for(int i=0; i<mVarLocations.size(); i++){
		*mVarLocations[i]=0.0;
	}
	for(int i=0; i<mInputLocations.size(); i++){
		*mInputLocations[i]=0.0;
	}
}

void iWQModel::setStateVariable(std::string name, double value)
{
	for(int i=0; i<mVarNames.size(); i++){
		if(mVarNames[i]==name){
			*mVarLocations[i]=value;
			break;
		}
	}
}

iWQStrings iWQModel::variableNames() const
{
	iWQStrings result;
	for(int i=0; i<mVarNames.size(); i++){
		if(!mShouldTakeDelta[i]){
			result.push_back(mVarNames[i]);
		}
	}
	return result;
}
//END DANGEROUS LOW-LEVEL FUNCTIONS

//############################################################################################################

#pragma mark CSTR channel transport models

iWQGenericChannelTransport::iWQGenericChannelTransport(std::string modelname) : iWQModel(modelname)
{
	mConcentrationVariable=NULL;
	
	//assume a unit cube
	channel_width = 1.0;
	channel_depth = 1.0;
	channel_length = 1.0;
	
	//make predefined stuff
	PAR(channel_width);
	PAR(channel_depth);
	PAR(channel_length);
	
	BFX(Fout);
	
	INP(Fin);
	INP(Fnew);
	INP(Qout);
}

//-------------------------------------------------------------------------------------------------------------

void iWQGenericChannelTransport::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	if(mConcentrationVariable){
		//actual concentration
		double C_val= *mConcentrationVariable;
		//concentration change		
		double V = channel_width * channel_depth * channel_length;
		*(D(mConcentrationVariable)) = (Fin + Fnew - Qout * C_val + R(x)) / V;
		//boundary flux
		*(F(Fout)) = Qout * C_val;
	}
}

//-------------------------------------------------------------------------------------------------------------

void iWQGenericChannelTransport::defineConcentrationVariable( double * var, const char * varname )
{
	if(!mConcentrationVariable){
		defineVariable(var, varname, false);
		mConcentrationVariable=var;
	}
}

//##################################################################################################

#pragma mark Initial values

iWQInitialValues::iWQInitialValues()
{
	mSharedParameterManager=NULL;
}

//--------------------------------------------------------------------------------------------------

bool iWQInitialValues::hasDefaultValueForVariable(std::string varname)
{
	//dynamic values
	if(mSharedParameterManager){	
		//try the initial value-specific dynamic value first
		std::string init_varname = varname + "@0";
		if(mSharedParameterManager->hasValueForParam(init_varname)){
			return true;
		}
		//then the dynamic but not specific
		if(mSharedParameterManager->hasValueForParam(varname)){
			return true;
		}
	}
	//local values
	return (mDefaultValues.find(varname)!=mDefaultValues.end());
}

//--------------------------------------------------------------------------------------------------

bool iWQInitialValues::hasValueForVariable(std::string varname, std::string modelid)
{
	bool result=false;
	//dynamic values
	if(mSharedParameterManager){	
		std::string init_varname = varname + "@0";
		//SCHEMA: variable@0[modelid]
		if(mSharedParameterManager->hasValueForParam(init_varname, modelid)){
			return true;
		}
		//SCHEMA: variable[modelid]
		if(mSharedParameterManager->hasValueForParam(varname, modelid)){
			return true;
		}
	}
	//local values
	if(mValues.find(modelid)!=mValues.end()){
		result=(mValues[modelid].find(varname)!=mValues[modelid].end());
	}
	else{
		//no entry for this modelid
		result=false; //hasDefaultValueForVariable(varname);
	}
	return result;
}

//--------------------------------------------------------------------------------------------------

double iWQInitialValues::defaultValueForVariable(std::string varname)
{
	if(mSharedParameterManager){	
		//try the initial value-specific dynamic value first
		std::string init_varname = varname + "@0";
		if(mSharedParameterManager->hasValueForParam(init_varname)){
			return mSharedParameterManager->valueForParam(init_varname);
		}
		//then the dynamic but not specific
		if(mSharedParameterManager->hasValueForParam(varname)){
			return mSharedParameterManager->valueForParam(varname);
		}
	}
	
	//then the static initial value at last (if any)
	return mDefaultValues[varname]; 
}

//--------------------------------------------------------------------------------------------------

double iWQInitialValues::valueForVariable(std::string varname, std::string modelid)
{
	if(mSharedParameterManager){	
		std::string init_varname = varname + "@0";
		//SCHEMA: variable@0[modelid]
		if(mSharedParameterManager->hasValueForParam(init_varname, modelid)){
			return mSharedParameterManager->valueForParam(init_varname, modelid);
		}
		//SCHEMA: variable[modelid]
		if(mSharedParameterManager->hasValueForParam(varname, modelid)){
			return mSharedParameterManager->valueForParam(varname, modelid);
		}
		//SCHEMA: variable@0
		if(mSharedParameterManager->hasValueForParam(init_varname)){
			return mSharedParameterManager->valueForParam(init_varname);
		}
		//SCHEMA: variable
		if(mSharedParameterManager->hasValueForParam(varname)){
			return mSharedParameterManager->valueForParam(varname);
		}
	}
	
	//local values at last
	if(hasValueForVariable(varname, modelid)){
		return mValues[modelid][varname];
	}
	else if(hasDefaultValueForVariable(varname)){
		return mDefaultValues[varname];
	}
	else{
		printf("[Error]: Requesting initial value for the unknown variable \"%s\".\n",varname.c_str());
		return 0.0;
	}
}

//--------------------------------------------------------------------------------------------------

void iWQInitialValues::setDefaultValueForVariable(double value, std::string varname)
{
	mDefaultValues[varname]=value;
}

//--------------------------------------------------------------------------------------------------

void iWQInitialValues::setValueForVariable(double value, std::string varname, std::string modelid)
{
	if(mValues.find(modelid)!=mValues.end()){
		mValues[modelid][varname]=value;
	}
	else{
		iWQKeyValues new_entry;
		new_entry[varname]=value;
		mValues[modelid]=new_entry;
	}
}

//--------------------------------------------------------------------------------------------------

iWQKeyValues iWQInitialValues::variablesForId(std::string modelid)
{
	if(mValues.find(modelid)!=mValues.end()){
		return mValues[modelid];
	}
	else{
		return mDefaultValues;
	}
}

//--------------------------------------------------------------------------------------------------

iWQKeyValues iWQInitialValues::defaultVariables()
{
	return mDefaultValues;
}

//--------------------------------------------------------------------------------------------------

iWQKeyValues iWQInitialValues::allValues()
{
	iWQKeyValues result=mDefaultValues;
	
	//add specific items
	std::map<std::string, iWQKeyValues>::iterator it1;
	for(it1=mValues.begin(); it1!=mValues.end(); ++it1){
		std::string modelid=it1->first;
		iWQKeyValues::iterator it2;
		iWQKeyValues vals=it1->second;
		for(it2=vals.begin(); it2!=vals.end(); ++it2){
			std::string key=it2->first;
			double value=it2->second;
			std::string merged_key = key + "[" + modelid + "]";
			result[merged_key]=value;
		}
	}
	return result;
}

//--------------------------------------------------------------------------------------------------

void iWQInitialValues::setParameterManager(iWQParameterManager * par)
{
	mSharedParameterManager = par;
}

//--------------------------------------------------------------------------------------------------

iWQParameterManager * iWQInitialValues::parameterManager()
{
	return mSharedParameterManager;
}
	
