/*
 *  evaluatormethod.cpp
 *  Various likelihood calculators
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */
 
//for the input-dependent bias method:
#include "Eigen/Dense"
#include "biasmatrices.h"

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include "evaluatormethod.h"
#include "model.h"
#include "datatable.h"

//-----------------------------------------------------------------------------------------------
#pragma mark Generic evaluator method functionality

iWQEvaluatorMethod::iWQEvaluatorMethod()
{ 
	mDataTable=0;
	mCommonParameters=0;
	initDefaultParams(); 
}

//-----------------------------------------------------------------------------------------------

bool iWQEvaluatorMethod::setParamValueFromMap(double * dest, std::string key, iWQSettingList * list, std::string flag)
{
	//convenience method to get a number from a settings map
	//printf("Searching for evaluation parameter %s.\n",key.c_str());
	bool result=false;
	if(hasDynamicParamValue(key, NULL, flag)){	//this checks if we have a value
		//prefer dynamic parameters over static ones
		setParamValueDynamically(dest, key, flag);	//this also updates it
		result=true;
		printf("Evaluation parameter %s for %s is updated dynamically.\n",key.c_str(), mComparisonLink.modelField().c_str());
	}
	
	//look for the static settings
	const char * p1;
	char * p2;
	if(key.length() && list && dest){
		std::string fullkey=key;
		bool hasflag=(flag.size()>0);
		if(hasflag){
			fullkey=key + "[" +flag + "]";	//if there is a valid flag
		}
		iWQSettingList::iterator it_full, it_plain;
		it_plain=list->find(key);
		if(hasflag){
			it_full=list->find(fullkey);
		}
		std::string valstr="";
		if(hasflag && it_full!=list->end()){
			//found the fully flagged parameter in the list
			valstr=it_full->second;
		}
		if(valstr.size()==0 && it_plain!=list->end()){
			valstr=it_plain->second;	//found it as a plain value
		}
		
		if(result && valstr.size()){
			//warn about conflicting definitions
			printf("[Warning]: Evaluation parameter %s was defined both as static and dynamic. Using dynamic mode.\n",fullkey.c_str());
		}
			
		//convert the number finally
		if(!result && valstr.length()){
			p1=valstr.c_str(); 
			double value=strtod(p1, &p2);
			if(p1!=p2){
				*dest=value;
				result=true;
				printf("Evaluation parameter %s for %s is static.\n",key.c_str(), mComparisonLink.modelField().c_str());
			}
		}
	}
		
	return result;
}

//-----------------------------------------------------------------------------------------------

bool iWQEvaluatorMethod::setParamValueFromMap(std::string * dest, std::string key, iWQSettingList * list, std::string flag)
{
	//convenience method to get a string from a settings map (no dymanical option, no flags effective)
	if(key.length() && list && dest){
		std::string fullkey=key;
		bool hasflag=(flag.size()>0);
		if(hasflag){
			fullkey=key + "[" +flag + "]";	//if there is a valid flag
		}
		iWQSettingList::iterator it_full, it_plain;
		it_plain=list->find(key);
		if(hasflag){
			it_full=list->find(fullkey);
		}
		std::string valstr="";
		if(hasflag && it_full!=list->end()){
			//found the fully flagged parameter in the list
			valstr=it_full->second;
		}
		if(valstr.size()==0 && it_plain!=list->end()){
			valstr=it_plain->second;	//found it as a plain value
		}
		
		//convert the number finally
		if(valstr.length()){
			*dest=valstr;
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------

bool iWQEvaluatorMethod::hasDynamicParamValue(std::string key, double * dest, std::string flag)		//responds with true if there is a properly named parameter and it updates the value too if dest is not NULL
{
	if(!mCommonParameters){
		printf("No parameter storage to browse from.\n");
		return false;
	}
	//specialized values first, general representation second (=prefer flagged if available)
	if(flag.size()==0){
		flag="EVAL";	//try this default flag for unspecified parameters
	}
	if(mCommonParameters->hasValueForParam(key, flag)){
		if(dest){
			*dest=mCommonParameters->valueForParam(key, flag);
		}
		return true;
	}
	else if(mCommonParameters->hasValueForParam(key)){
		if(dest){
			*dest=mCommonParameters->valueForParam(key);
		}
		return true;
	}
	return false;	//no matching parameter found in dynamic storage
}

//-----------------------------------------------------------------------------------------------

bool iWQEvaluatorMethod::setParamValueDynamically(double * dest, std::string key, std::string flag)
{
	if(dest && hasDynamicParamValue(key,dest,flag)){
		std::string fullkey=key;
		if(flag.size()){
			//code flag into key in an easily separable manner
			fullkey=key+" "+flag;	//space is not allowed in XML names, so it's good for separation
		}
		mDynamicParams[fullkey]=dest;
		return true;
	}
	else{
		return false;
	}
}

//-----------------------------------------------------------------------------------------------

void iWQEvaluatorMethod::updateDynamicParams()	//will be called by the evaluator from outside on parameter update	
{
	std::map<std::string, double *>::iterator it;
	for(it=mDynamicParams.begin(); it!=mDynamicParams.end(); ++it){
		std::string fullkey = it->first;
		std::string key=fullkey;
		std::string flag="";
		//split key if necessary
		size_t pos=fullkey.find(" ");
		if(pos!=std::string::npos){
			key=fullkey.substr(0,pos);
			flag=fullkey.substr(pos+1);
		}
		double * dest = it->second;
		if(dest){
			hasDynamicParamValue(key, dest, flag);
		}
	}
}

//-----------------------------------------------------------------------------------------------

double iWQEvaluatorMethod::evaluate()
{
	if(!mDataTable){
		return DBL_MAX;	//bad by default
	}

	int start=0;
	int end=mDataTable->numRows();
	return evaluate(start, end);
}

//===============================================================================================

#pragma mark Class factory function

iWQEvaluatorMethod * createEvalMethod(std::string methodName)
{
	iWQEvaluatorMethod * evalMethod=NULL;
	
	if(methodName.compare("Nash-Sutcliffe")==0 || methodName.compare("NS")==0 || methodName.compare("NSBoxCox")==0){
		//alloc a simple NS method
		evalMethod=new iWQNSBoxCoxEvaluation;
	}
	else if(methodName.compare("Normal Error")==0 || methodName.compare("LogLikeliNormal")==0 || methodName.compare("Normal")==0){
		//alloc a normally distributed error likelihood method
		evalMethod=new iWQNormalLikelihoodEvaluation;
	}
	else if(methodName.compare("Heteroscedastic Normal Error")==0 || methodName.compare("LogLikeliHetNormal")==0 || methodName.compare("HetNormal")==0){
		//alloc a normally distributed error likelihood method
		evalMethod=new iWQHeteroscedasticNormalLikelihoodEvaluation;
	}
	else if(methodName.compare("Quantile Normal Error")==0 || methodName.compare("LogLikeliQuantileNormal")==0 || methodName.compare("QuantileNormal")==0 || methodName.compare("QuantNormal")==0 || methodName.compare("QNormal")==0){
		//alloc a normally distributed quantiles error likelihood method
		evalMethod=new iWQQuantileNormalLikelihoodEvaluation;
	}
	else if(methodName.compare("Quantile Error")==0 || methodName.compare("LogLikeliQuantile")==0 || methodName.compare("Quantile")==0 || methodName.compare("Quant")==0 || methodName.compare("Q")==0){
		//alloc a normally distributed quantiles error likelihood method
		evalMethod=new iWQQuantileLikelihoodEvaluation;
	}
	else if(methodName.compare("Input-dependent Bias")==0 || methodName.compare("LogLikeliIDAR")==0 || methodName.compare("IDAR")==0){
		//alloc an IDAR likelihood method
		evalMethod=new iWQIDARLikelihoodEvaluation;
	}
	else if(methodName.compare("Input-dependent Bias and Normal Error")==0 || methodName.compare("LogLikeliBIAS")==0 || methodName.compare("BIAS")==0){
		//alloc a BIAS-IDAR likelihood method
		evalMethod=new iWQBiasIDARLikelihoodEvaluation;
	}
	else if(methodName.compare("AR1 with SEP innovations")==0 || methodName.compare("LogLikeliARSEP")==0 || methodName.compare("ARSEP")==0){
		//alloc a BIAS-IDAR likelihood method
		evalMethod=new iWQARSEPLikelihoodEvaluation;
	}
			
	if(evalMethod){
		evalMethod->initDefaultParams();
	}
	return evalMethod;
}

//===============================================================================================

#pragma mark Nash-Sutcliffe index with Box-Cox transformation

double iWQNSBoxCoxEvaluation::evaluate(int startindex, int endindex)
{
	//Nash-Sutcliffe statistics on the Box-Cox transformed values
	//returns NaN or INF if the transformation fails for any value
	double sum=0.0;
	double sumsqdeviation=0.0;
	double sumsqmodeldeviation=0.0;
	int count=0;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			sum+=mComparisonLink.measurement();	//<----
			count++;
		}
	}
	
	double average=sum/(double)count;
		
	//deviations from the averages (Nash-Sutcliffe statistics)
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=boxcox_transform(lambda_1, lambda_2, mComparisonLink.measurement(), NULL);		//<----
			double model=boxcox_transform(lambda_1, lambda_2, mComparisonLink.model(), NULL);			//<----
			sumsqdeviation+=(meas-average)*(meas-average);
			sumsqmodeldeviation+=(meas-model)*(meas-model);
		}
	}
	
	double NS=(sumsqdeviation!=0.0)?sumsqmodeldeviation/sumsqdeviation:0.0;
	
	return NS;	
}

//--------------------------------------------------------------------------------------------------

void iWQNSBoxCoxEvaluation::setParams(iWQSettingList list)
{
	//in-place settings
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&lambda_1,"lambda_1",&list,varname);
	setParamValueFromMap(&lambda_2,"lambda_2",&list,varname);
} 

//--------------------------------------------------------------------------------------------------

void iWQNSBoxCoxEvaluation::initDefaultParams()
{ 
	lambda_1=1.0; 
	lambda_2=0.0; 
}
//--------------------------------------------------------------------------------------------------

#pragma mark i.i.d. normal error model 

void iWQNormalLikelihoodEvaluation::initDefaultParams()
{
	//defaults to standard normally distributed errors
	sigma=1.0;
	dist.setMean(0);
	lambda_1=1.0;
	lambda_2=0.0;
	LOQ=-DBL_MAX;
	sumlogy = -DBL_MAX;
}	

void iWQNormalLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma,"sigma",&list,varname);
	setParamValueFromMap(&lambda_1,"lambda_1",&list,varname);
	setParamValueFromMap(&lambda_2,"lambda_2",&list,varname);
	setParamValueFromMap(&LOQ,"LOQ",&list,varname);
}

double iWQNormalLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with normal error model
	double loglikeli=0.0;
	double result=0.0;
	
	dist.setStdev(sigma);	//update before each evaluation
	
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int j=startindex; j<endindex; j++){
			mDataTable->setRow(j);
			if(mComparisonLink.numeric()){
				double meas_raw = mComparisonLink.measurement();
				if(meas_raw<=LOQ){
					meas_raw = 0.5 * LOQ;
				}
				if(meas_raw + lambda_2 > 0.0){
					sumlogy += log(meas_raw + lambda_2);
				}
				else{
					printf("[Warning]: Measurement (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), meas_raw, j);
					sumlogy=0.0;
					break;
				}
			}
		}
	}	
		
	loglikeli += (lambda_1 - 1.0) * sumlogy;
		
	//log likelihood of deviations
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas_raw = mComparisonLink.measurement();
			double model_raw = mComparisonLink.model();
			if(meas_raw>LOQ){ //non-LOQ measurements
				double meas=boxcox_transform(lambda_1, lambda_2, meas_raw, NULL);
				double model=boxcox_transform(lambda_1, lambda_2, model_raw, NULL);
				double newres=dist.logLikeli(model-meas);
				loglikeli+=newres;
				if(std::isnan(newres) || std::isinf(newres) || newres== DBL_MAX || newres== -DBL_MAX){
					printf("[Warning]: Log likelihood of point %d in %s (measured=%lf, modelled=%lf) is %lf\n",j, modelFieldName().c_str(), meas_raw, model_raw, newres);
					printf("           BC(measured)=%lf, BC(modelled)=%lf\n", meas, model);
					printf("           lambda_1=%lf, lambda_2=%lf\n", lambda_1, lambda_2);
				}
			}
			else{
				//get the cumulative likelihood that meas_raw is below LOQ
				double meas=boxcox_transform(lambda_1, lambda_2, LOQ, NULL);
				double model=boxcox_transform(lambda_1, lambda_2, model_raw, NULL);
				double xi = (meas-model) / sigma; //standardized variable
				double pxi = lpnorm(xi);
				double pxinull = 0.0;
				double newres = pxi; //log(pxi - pxinull);
				if(std::isnan(newres) || std::isinf(newres) || newres== DBL_MAX || newres== -DBL_MAX){
					printf("[Warning]: Log likelihood of point %d in %s (measured=%lf (<=LOQ), modelled=%lf) is %lf\n",j, modelFieldName().c_str(), meas_raw, model_raw, newres);
					printf("           BC(measured)=%lf, BC(modelled)=%lf\n", meas, model);
					printf("           lambda_1=%lf, lambda_2=%lf\n", lambda_1, lambda_2);
					printf("           sigma=%lf, xi=%lf, p(xi)=%lf, p(xi0)=%lf\n", sigma, xi, pxi, pxinull);
				}
				loglikeli+=newres;
			}
		}
	}
	return -loglikeli;	//to make it reversed for minimization
}

std::vector<std::string> iWQNormalLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	
	result.push_back(std::string("Y_")+varname);
	result.push_back(std::string("YE_")+varname);
	result.push_back(std::string("Ytr_")+varname);
	result.push_back(std::string("YEtr_")+varname);
	
	return result;
}

void iWQNormalLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	std::vector<double> Ys;
	std::vector<double> Ytrs;
	std::vector<double> YEs;
	std::vector<double> YEtrs;
	
	dist.setStdev(sigma);
	
	mDataTable->rewind();
	while(mDataTable->stepRow()!=-1){
		double model=mComparisonLink.model();
		double modeltr = boxcox_transform(lambda_1, lambda_2, model, NULL);
		Ys.push_back(model);
		Ytrs.push_back(modeltr);
		if(mComparisonLink.numeric()){
			//past
			double meas=mComparisonLink.measurement();
			double meastr = boxcox_transform(lambda_1, lambda_2, meas, NULL);
			YEs.push_back(meas);
			YEtrs.push_back(meastr);
		}
		else{
			dist.setStdev(sigma);
			double Etr = dist.generate();
			double YEtr = modeltr + Etr;
			double YE = boxcox_retransform(lambda_1, lambda_2, YEtr, NULL);
			YEtrs.push_back(YEtr);
			YEs.push_back(YE);
		}
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("Ytr_"+varname)=Ytrs;
	storage->operator[]("YE_"+varname)=YEs;
	storage->operator[]("YEtr_"+varname)=YEtrs;
}

//----------------------------------------------------------------------------------------

#pragma mark Heteroscedastic normal error model with transformation

void iWQHeteroscedasticNormalLikelihoodEvaluation::initDefaultParams()
{
	iWQNormalLikelihoodEvaluation::initDefaultParams();
	inputfieldname="";
	inputptr=NULL;
	k_input=1.0;
}	
	
void iWQHeteroscedasticNormalLikelihoodEvaluation::setParams(iWQSettingList list)
{
	iWQNormalLikelihoodEvaluation::setParams(list);
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&inputfieldname,"driver",&list,varname);	
	setParamValueFromMap(&k_input,"k_input",&list,varname);
}

double iWQHeteroscedasticNormalLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with normal error model
	double loglikeli=0.0;
	double result=0.0;
	
	dist.setStdev(sigma);	//update before each evaluation
	inputptr=mDataTable->portForColumn(inputfieldname);
		
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int j=startindex; j<endindex; j++){
			mDataTable->setRow(j);
			if(mComparisonLink.numeric()){
				double meas_raw = mComparisonLink.measurement();
				if(meas_raw<=LOQ){
					meas_raw = 0.5 * LOQ;
				}
				if(meas_raw + lambda_2 > 0.0){
					sumlogy += log(meas_raw + lambda_2);
				}
				else{
					printf("[Warning]: Measurement (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), meas_raw, j);
					sumlogy=0.0;
					break;
				}
			}
		}
	}	
		
	loglikeli += (lambda_1 - 1.0) * sumlogy;
		
	//log likelihood of deviations
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas_raw = mComparisonLink.measurement();
			double model_raw = mComparisonLink.model();
			double act_scaling = (inputptr && *inputptr!=DBL_MAX && *inputptr>0.0 && k_input>0.0) ? *inputptr/k_input : 1.0;
			if(meas_raw>LOQ){ //non-LOQ measurements
				double meas=boxcox_transform(lambda_1, lambda_2, meas_raw, NULL);
				double model=boxcox_transform(lambda_1, lambda_2, model_raw, NULL);
				double newres=dist.logLikeli((model-meas)/act_scaling);
				loglikeli+=newres;
				if(std::isnan(newres) || std::isinf(newres) || newres== DBL_MAX || newres== -DBL_MAX){
					printf("[Warning]: Log likelihood of point %d in %s (measured=%lf, modelled=%lf) is %lf\n",j, modelFieldName().c_str(), meas_raw, model_raw, newres);
					printf("           BC(measured)=%lf, BC(modelled)=%lf\n", meas, model);
					printf("           lambda_1=%lf, lambda_2=%lf\n", lambda_1, lambda_2);
				}
			}
			else{
				//get the cumulative likelihood that meas_raw is below LOQ
				double meas=boxcox_transform(lambda_1, lambda_2, LOQ, NULL);
				double model=boxcox_transform(lambda_1, lambda_2, model_raw, NULL);
				double xi = (meas-model) / (sigma * act_scaling); //standardized variable
				double pxi = lpnorm(xi);
				double pxinull = 0.0;
				double newres = pxi; //log(pxi - pxinull);
				if(std::isnan(newres) || std::isinf(newres) || newres== DBL_MAX || newres== -DBL_MAX){
					printf("[Warning]: Log likelihood of point %d in %s (measured=%lf (<=LOQ), modelled=%lf) is %lf\n",j, modelFieldName().c_str(), meas_raw, model_raw, newres);
					printf("           BC(measured)=%lf, BC(modelled)=%lf\n", meas, model);
					printf("           lambda_1=%lf, lambda_2=%lf\n", lambda_1, lambda_2);
					printf("           sigma=%lf, xi=%lf, p(xi)=%lf, p(xi0)=%lf\n", sigma, xi, pxi, pxinull);
				}
				loglikeli+=newres;
			}
		}
	}
	return -loglikeli;	//to make it reversed for minimization
}

void iWQHeteroscedasticNormalLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	std::vector<double> Ys;
	std::vector<double> Ytrs;
	std::vector<double> YEs;
	std::vector<double> YEtrs;
	
	dist.setStdev(sigma);
	inputptr=mDataTable->portForColumn(inputfieldname);
	
	mDataTable->rewind();
	while(mDataTable->stepRow()!=-1){
		double act_scaling = (inputptr && *inputptr!=DBL_MAX && *inputptr>0.0 && k_input>0.0) ? *inputptr/k_input : 1.0;
		double model=mComparisonLink.model();
		double modeltr = boxcox_transform(lambda_1, lambda_2, model, NULL);
		Ys.push_back(model);
		Ytrs.push_back(modeltr);
		if(mComparisonLink.numeric()){
			//past
			double meas=mComparisonLink.measurement();
			double meastr = boxcox_transform(lambda_1, lambda_2, meas, NULL);
			YEs.push_back(meas);
			YEtrs.push_back(meastr);
		}
		else{
			dist.setStdev(sigma);
			double Etr = dist.generate();
			double YEtr = modeltr + Etr * act_scaling;
			double YE = boxcox_retransform(lambda_1, lambda_2, YEtr, NULL);
			YEtrs.push_back(YEtr);
			YEs.push_back(YE);
		}
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("Ytr_"+varname)=Ytrs;
	storage->operator[]("YE_"+varname)=YEs;
	storage->operator[]("YEtr_"+varname)=YEtrs;
}

//--------------------------------------------------------------------------------------------------

#pragma mark Normal error model on the series quantiles

void iWQQuantileNormalLikelihoodEvaluation::populateProbs()
{
	probs.clear();
	int n_ps_half = (int)(0.5/quantspacing);
	for(int i=-n_ps_half; i<=n_ps_half; i++){
		double p = 0.5+(double)i * quantspacing;
		if(p>0 && p<1){
			probs.push_back(p);
		}
	}
}

void iWQQuantileNormalLikelihoodEvaluation::initDefaultParams()
{
	//defaults to standard normally distributed errors
	sigma=1.0;
	dist.setMean(0);
	lambda_1=1.0;			//no transfromation
	lambda_2=0.0;
	quantspacing = 0.475; 	//3 default quantiles: 0.025, 0.5, 0.975
	LOQ=-DBL_MAX;
	sumlogy=-DBL_MAX;
	populateProbs();
}	

void iWQQuantileNormalLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma,"sigma",&list,varname);
	setParamValueFromMap(&lambda_1,"lambda_1",&list,varname);
	setParamValueFromMap(&lambda_2,"lambda_2",&list,varname);
	setParamValueFromMap(&LOQ,"LOQ",&list,varname);
	setParamValueFromMap(&quantspacing,"quantspacing",&list,varname);
	if(quantspacing<=0.0){
		quantspacing=0.475;	//revert to default for 0 or negative spacing
	}
	//update the target probabilities
	populateProbs();
}

double iWQQuantileNormalLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with normal error model
	double loglikeli=0.0;
	
	//prepare the lists of measurements and models
	std::vector<double> measured;
	std::vector<double> modelled;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=boxcox_transform(lambda_1, lambda_2, mComparisonLink.measurement(), NULL);
			double model=boxcox_transform(lambda_1, lambda_2, mComparisonLink.model(), NULL);
			measured.push_back(meas);
			modelled.push_back(model);
		}
	}
	
	std::sort(measured.begin(), measured.end());
	std::sort(modelled.begin(), modelled.end());
	
	//log likelihood of quantile deviations
	//prepare quantiles
	std::vector<double> q_hat;
	std::vector<double> q;
	for(int i=0; i<probs.size(); i++){
		q_hat.push_back(quantile(measured, probs[i], 7, true));
		q.push_back(quantile(modelled, probs[i], 7, true));
	}
	//check LOQ position
	int startpos = 0;
	if(LOQ!=-DBL_MAX){
		double LOQtr = boxcox_transform(lambda_1, lambda_2, LOQ, NULL);
		for(int i=0; i<probs.size(); i++){
			if(q_hat[i]>LOQtr){
				startpos=i>0?i-1:0;
				break;
			}
		}
	}
	
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int i=startpos; i<q_hat.size(); i++){
			if(q_hat[i] + lambda_2 > 0.0){
				sumlogy += log(q_hat[i] + lambda_2);
			}
			else{
				printf("[Warning]: Measurement quantile (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), q_hat[i], i);
				sumlogy=0.0;
				break;
			}
		}
	}
	
	loglikeli += (lambda_1 - 1.0) * sumlogy;
	
	//now assess the likelihood from startpos
	for(int i=startpos; i<probs.size(); i++){
		double p_uncond_upper = 1.0;
		double p_uncond_lower = 0.0;
		if(i>startpos){
			p_uncond_lower = pnorm((q_hat[i-1] - q[i])/sigma); //0.5 * (1 + erf((q_hat[i-1] - q[i])/(1.414213562 * sigma)));	//P(q_hat[i] < q_hat[i-1])
		}
		if(i<probs.size()-1){
			p_uncond_upper = pnorm((q_hat[i+1] - q[i])/sigma); //0.5 * (1 + erf((q_hat[i+1] - q[i])/(1.414213562 * sigma)));	//P(q_hat[i] > q_hat[i+1])
		}
		double p_cond = p_uncond_upper - p_uncond_lower; //1.0 - ;	//Normal CDF for q1hat with mu=q2 and sd=sigma
		dist.setMean(q[i]);
		dist.setStdev(sigma);
		double ll_uncond = dist.logLikeli(q_hat[i]);
		if(p_cond>0.0 && !std::isinf(ll_uncond)){
			loglikeli += ll_uncond - log(p_cond);
		}
		else{
			return 0.99*DBL_MAX;	//ll_uncond seems to be 0.0, so return _almost_ INF (this is not the likelihood=0 case but we don't have enough numerical accuracy to calculate the likelihood)
		}
		//}
		//q1hat = q_meas;
	}
	return -loglikeli;	//to make it reversed for minimization
}

std::vector<std::string> iWQQuantileNormalLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	
	result.push_back(std::string("Q_")+varname);	//quantiles in normal space
	result.push_back(std::string("Qtr_")+varname);	//quantiles in transformed space
	result.push_back(std::string("QE_")+varname);	//predictive Q+E in normal space
	result.push_back(std::string("QEtr_")+varname);	//predictive Q+E in transformed space
	result.push_back(std::string("Y_")+varname);	//series in normal space
	result.push_back(std::string("Ytr_")+varname);	//series in transformed space
		
	return result;
}

void iWQQuantileNormalLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	std::vector<double> Ys;
	std::vector<double> Ytrs;
	std::vector<double> Qs;
	std::vector<double> Qtrs;
	std::vector<double> QEs;
	std::vector<double> QEtrs;
		
	mDataTable->rewind();
	while(mDataTable->stepRow()!=-1){
		double model=mComparisonLink.model();
		double modeltr = boxcox_transform(lambda_1, lambda_2, model, NULL);
		Ys.push_back(model);
		Ytrs.push_back(modeltr);
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("Ytr_"+varname)=Ytrs;
		
	//make quantiles
	std::vector<double> modelled = Ys;
	std::vector<double> modelled_tr = Ytrs;
	std::sort(modelled.begin(), modelled.end());
	std::sort(modelled_tr.begin(), modelled_tr.end());
	
	for(int i=0; i<probs.size(); i++){
		double p = probs[i];
		Qs.push_back(quantile(modelled, p, 7, true));
		Qtrs.push_back(quantile(modelled_tr, p, 7, true));
	}
	
	storage->operator[]("Q_"+varname)=Qs;
	storage->operator[]("Qtr_"+varname)=Qtrs;
	
	//do predictive Q+E with Gibbs sampling(1 realisation)
	int ngibbs = 500;	//hardcoded gibbs sample length
	
	QEtrs = Qtrs;
	int nqs = Qtrs.size();
	
	for(int j=0; j<=ngibbs; j++){
		for(int i=0; i<nqs; i++){
			//draw a QE from the conditional likelihood
			double * cond_low = NULL;
			double * cond_high = NULL;
			if(i>0){
				cond_low = &QEtrs[i-1]; 
			}
			if(i<nqs-1){ 
				cond_high = &QEtrs[i+1]; 
			}
			QEtrs[i] = rtnorm(Qtrs[i], sigma, cond_low, cond_high);
		}
	}
	//retransform QEtrs into normal space
	QEs=Qs;	//lazy way of allocation
	for(int i=0; i<nqs; i++){
		QEs[i] = boxcox_retransform(lambda_1, lambda_2, QEtrs[i], NULL);
	}
	
	storage->operator[]("QE_"+varname)=QEs;
	storage->operator[]("QEtr_"+varname)=QEtrs;
	
}
 
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------

#pragma mark Quantile error model version 2

void iWQQuantileLikelihoodEvaluation::populateProbs()
{
	probs.clear();
	int n_ps_half = (int)(0.5/quantspacing);
	for(int i=-n_ps_half; i<=n_ps_half; i++){
		double p = 0.5+(double)i * quantspacing;
		if(p>0 && p<1){
			probs.push_back(p);
		}
	}
}

void iWQQuantileLikelihoodEvaluation::initDefaultParams()
{
	//defaults to standard normally distributed errors
	sigma=1.0;
	dist.setMean(0);
	lambda_1=1.0;			//no transformation
	lambda_2=0.0;
	quantspacing = 0.475; 	//3 default quantiles: 0.025, 0.5, 0.975
	LOQ=-DBL_MAX;
	populateProbs();
}	

void iWQQuantileLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma,"sigma",&list,varname);
	setParamValueFromMap(&LOQ,"LOQ",&list,varname);
	setParamValueFromMap(&quantspacing,"quantspacing",&list,varname);
	if(quantspacing<=0.0){
		quantspacing=0.475;	//revert to default for 0 or negative spacing
	}
	//update the target probabilities
	populateProbs();
}

double iWQQuantileLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with normal error model
	double loglikeli=0.0;
	
	//prepare the lists of measurements and models
	std::vector<double> measured;
	std::vector<double> modelled;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=mComparisonLink.measurement();
			double model=mComparisonLink.model();
			measured.push_back(meas);
			modelled.push_back(model);
		}
	}
	
	std::sort(measured.begin(), measured.end());
	std::sort(modelled.begin(), modelled.end());
	
	//log likelihood of quantile deviations
	//prepare quantiles
	std::vector<double> q_hat;
	std::vector<double> q;
	for(int i=0; i<probs.size(); i++){
		q_hat.push_back(quantile(measured, probs[i], 7, true));
		q.push_back(quantile(modelled, probs[i], 7, true));
	}
	//check LOQ position
	int startpos = 0;
	if(LOQ!=-DBL_MAX){
		double LOQtr = boxcox_transform(lambda_1, lambda_2, LOQ, NULL);
		for(int i=0; i<probs.size(); i++){
			if(q_hat[i]>LOQtr){
				startpos=i>0?i-1:0;
				break;
			}
		}
	}
	
	//prepare sigma vector
	std::vector<double> densities = density(modelled, q);
	
	if(densities.size()!=q.size()){
		return DBL_MAX;
	}
	
	//now assess the likelihood from startpos
	for(int i=startpos; i<probs.size(); i++){
		dist.setMean(q[i]);
		double densi = densities[i];
		double idensi2 = 1.0 / (densi * densi);
		dist.setStdev(sqrt(sigma * probs[i] * (1.0 - probs[i]) * idensi2));
		loglikeli += dist.logLikeli(q_hat[i]);
	}
	return -loglikeli;	//to make it reversed for minimization
}

std::vector<std::string> iWQQuantileLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	
	result.push_back(std::string("Q_")+varname);	//quantiles in normal space
	result.push_back(std::string("QE_")+varname);	//predictive Q+E in normal space
	result.push_back(std::string("Y_")+varname);	//series in normal space
		
	return result;
}

void iWQQuantileLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	std::vector<double> Ys;
	std::vector<double> Qs;
	std::vector<double> QEs;
		
	mDataTable->rewind();
	while(mDataTable->stepRow()!=-1){
		double model=mComparisonLink.model();
		double modeltr = model;
		Ys.push_back(model);
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
		
	//make quantiles
	std::vector<double> modelled = Ys;
	std::sort(modelled.begin(), modelled.end());
	
	for(int i=0; i<probs.size(); i++){
		double p = probs[i];
		Qs.push_back(quantile(modelled, p, 7, true));
	}
	
	storage->operator[]("Q_"+varname)=Qs;
	
	//do predictive Q+E
	QEs = Qs;
	int nqs = Qs.size();
	
	//prepare sigma vector
	std::vector<double> densities = density(modelled, Qs);
	
	for(int i=0; i<nqs; i++){
		//draw a QE
		double densi = densities[i];
		double idensi2 = 1.0 / (densi * densi);
		double sigmaq = sigma * probs[i] * (1.0 - probs[i]) * idensi2;
		dist.setMean(Qs[i]);
		dist.setStdev(sigmaq);
		QEs[i] = dist.generate();
	}
	
	storage->operator[]("QE_"+varname)=QEs;
}

//--------------------------------------------------------------------------------------------------

#pragma mark Input-dependent AutoRegressive (IDAR) error model

void iWQIDARLikelihoodEvaluation::initDefaultParams()
{
	//inits the built-in distribution
	dist.setMean(0.0);
	dist.setStdev(1.0);
	sigma_b2=1.0;
	beta=20.0;	//practically no correlation
	kappa=0.0;
		
	lambda_1 = 1.0;	//defaults to no transform
	lambda_2 = 0.0;
	
	inputfieldname="";
	inputptr=NULL;
	sumlogy = -DBL_MAX;
}

void iWQIDARLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma_b2,"sigma_b2",&list, varname);
	setParamValueFromMap(&beta,"beta",&list, varname);
	setParamValueFromMap(&kappa,"kappa",&list, varname);
	setParamValueFromMap(&lambda_1,"lambda_1",&list, varname);	//transformation parameters
	setParamValueFromMap(&lambda_2,"lambda_2",&list, varname);
	
	setParamValueFromMap(&inputfieldname,"driver",&list, varname);
}

double iWQIDARLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with IDAR error model
	double loglikeli=0.0;
		
	inputptr=mDataTable->portForColumn(inputfieldname);
	
	if(!inputptr){
		return loglikeli;
	}
	
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int j=startindex; j<endindex; j++){
			mDataTable->setRow(j);
			if(mComparisonLink.numeric()){
				double meas_raw = mComparisonLink.measurement();
				if(meas_raw + lambda_2 > 0.0){
					sumlogy += log(meas_raw + lambda_2);
				}
				else{
					printf("[Warning]: Measurement (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), meas_raw, j);
					sumlogy=0.0;
					break;
				}
			}
		}
	}	
		
	loglikeli += (lambda_1 - 1.0) * sumlogy;
	
	//log likelihood of deviations
	double rho = exp(-beta);
	
	double prev_bias=0.0;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=mComparisonLink.measurement();
			double model=mComparisonLink.model();
			double meastr=boxcox_transform(lambda_1, lambda_2, meas, NULL);
			double modeltr=boxcox_transform(lambda_1, lambda_2, model, NULL);
			double act_bias=modeltr-meastr;	//was model-meas
			double input=*inputptr;
			
			//reformulated
			double condstdev = sqrt(jumpVarianceOfB(sigma_b2, beta, kappa, 0.0, input));
			double condmean = rho * prev_bias;
			dist.setMean(condmean);
			dist.setStdev(condstdev);
			loglikeli+=dist.logLikeli(act_bias);
			
			prev_bias=act_bias;
		}
	}
	return -loglikeli;	//to make it reversed for minimization
}

std::vector<std::string> iWQIDARLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	result.push_back("Y_"+varname);
	result.push_back("YB_"+varname);
	result.push_back("Ytr_"+varname);
	result.push_back("YBtr_"+varname);
	result.push_back("I_"+varname);
	return result;
}

void iWQIDARLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	std::vector<double> Ys;
	std::vector<double> YBs;
	std::vector<double> Ytrs;
	std::vector<double> YBtrs;
	std::vector<double> Is;
	
	inputptr=mDataTable->portForColumn(inputfieldname);
	
	double rho = exp(-beta);
		
	mDataTable->rewind();
	double prev_bias=0.0;
	while(mDataTable->stepRow()!=-1){
		double model=mComparisonLink.model();
		double modeltr=boxcox_transform(lambda_1, lambda_2, model, NULL);
		Ys.push_back(model);
		Ytrs.push_back(modeltr);
		if(mComparisonLink.numeric()){
			//past
			double meas=mComparisonLink.measurement();
			double meastr=boxcox_transform(lambda_1, lambda_2, meas, NULL);
			double act_bias=modeltr-meastr;
			YBs.push_back(meas);
			YBtrs.push_back(meastr);
			
			//random increment of the past
			double input=*inputptr;
			double condstdev = sqrt(jumpVarianceOfB(sigma_b2, beta, kappa, 0.0, input));
			double condmean = rho * prev_bias;
			double jump = (act_bias - condmean)/(condstdev!=0.0?condstdev:1.0);
			Is.push_back(jump);
			
			prev_bias=act_bias;
		}
		else{
			//prediction stage
			double input=*inputptr;
			
			double condstdev = sqrt(jumpVarianceOfB(sigma_b2, beta, kappa, 0.0, input));
			double condmean = rho * prev_bias;
			dist.setMean(condmean);
			dist.setStdev(condstdev);
			double val=dist.generate();	//bias with transformation
						
			YBs.push_back(boxcox_retransform(lambda_1, lambda_2, modeltr-val, NULL));
			YBtrs.push_back(modeltr-val);
						
			double jump = (val - condmean)/(condstdev!=0.0?condstdev:1.0);
			Is.push_back(jump);
			
			prev_bias=val;
		}
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("YB_"+varname)=YBs;
	storage->operator[]("Ytr_"+varname)=Ytrs;
	storage->operator[]("YBtr_"+varname)=YBtrs;
	storage->operator[]("I_"+varname)=Is;
}

//------------------------------------------------------------------

#pragma mark IDAR model bias and independent measurement error

void iWQBiasIDARLikelihoodEvaluation::initDefaultParams()
{
	//inits the built-in distribution
	dist.setMean(0.0);
	dist.setStdev(1.0);
	//default: unit bias (with tcorr=1) and unit noise variance
	sigma_b2=1.0;
	beta=20.0;
	kappa=0.0;
	sigma_e2=1.0;
	inputfieldname="";
	inputptr=NULL;
	pi = 0.0;
	kappa_e = 0.0;
	
	lambda_1=1.0;
	lambda_2=0.0;
	sumlogy = -DBL_MAX;
	
	maxkernelsize = 10;		//must be even
}

void iWQBiasIDARLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma_b2,"sigma_b2",&list,varname);
	setParamValueFromMap(&beta,"beta",&list,varname);
	setParamValueFromMap(&sigma_e2,"sigma_e2",&list,varname);
	setParamValueFromMap(&kappa,"kappa",&list,varname);
	setParamValueFromMap(&pi,"pi",&list,varname);
	setParamValueFromMap(&kappa_e,"kappa_e",&list,varname);
	
	setParamValueFromMap(&inputfieldname,"driver",&list,varname);
	
	setParamValueFromMap(&lambda_1,"lambda_1",&list,varname);	//transformation parameters
	setParamValueFromMap(&lambda_2,"lambda_2",&list,varname);
	
	setParamValueFromMap(&maxkernelsize,"max_kernel_size",&list,varname);
}

double iWQBiasIDARLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//pre-filter parameters
	double minbeta=1E-3;
	double maxbeta=10.0;
	double minsigma=1E-8;
	double minkappa=0.0;
	
	//std::cout<<"Evaluate START"<<std::endl;
	
	int penalty=0;
	if(beta<minbeta || beta>maxbeta){
		penalty++;
	}
	if(sigma_e2<minsigma || sigma_b2<minsigma){
		penalty++;
	}
	if(kappa<minkappa){
		penalty++;
	}
	if(kappa_e<minkappa){
		penalty++;
	}
	
	if(penalty>0){
		return DBL_MAX;
	}
	
	//make residual series
	inputptr=mDataTable->portForColumn(inputfieldname);
	std::vector<double> yL_yLM;
	std::vector<double> inputs;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=boxcox_transform(lambda_1, lambda_2, mComparisonLink.measurement(), NULL);
			double model=boxcox_transform(lambda_1, lambda_2, mComparisonLink.model(), NULL);
			yL_yLM.push_back(meas-model);
			inputs.push_back(*inputptr);
		}
		else{
			break;
		}
	}
	int dim=yL_yLM.size();
	
	//calculate likelihood
	double result=0.0;
	double loglikeli=0.0;
	double pipart = log(1.0 / sqrt(2.0 * M_PI)); 
	
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int j=startindex; j<endindex; j++){
			mDataTable->setRow(j);
			if(mComparisonLink.numeric()){
				double meas_raw = mComparisonLink.measurement();
				if(meas_raw + lambda_2 > 0.0){
					sumlogy += log(meas_raw + lambda_2);
				}
				else{
					printf("[Warning]: Measurement (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), meas_raw, j);
					sumlogy=0.0;
					break;
				}
			}
		}
	}	
		
	loglikeli += (lambda_1 - 1.0) * sumlogy;
	
	int md=(int)maxkernelsize;	//fixed kernel size, was 10
	if(md%1){
		md++;
	}
	if(md<4){
		md=4;
	}
	if(md<dim){
		//kernel solution
		double det, det1;
		int md1=md+1;
		std::vector<double> inp (md);
		std::vector<double> inp1 (md1);
		
		Eigen::MatrixXd kernelinv; 
		Eigen::MatrixXd kernelinv1;
		
		//serially evaluate the likelihood
		//process the residuals by moving the kernels around them
		Eigen::VectorXd window1 (md1);		//md+1 elements
		Eigen::VectorXd window (md);		//md elements
		
		double exppart, exppart1;
		double lik, lik1;
		
		//process them according to conditional probability
		for(int j=0; j<=dim-md1; j++){
			//take md1 elements from the output
			for(int i=0; i<md1; i++){
				window1[i]=yL_yLM[j+i];
				inp1[i]=inputs[j+i];
			}
					
			//create the outer kernel
			kernelinv1=makeCovarMatrix(inp1, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e, &det1);
			exppart1=(window1.transpose() * kernelinv1).dot(window1);
			lik1 = md1 * pipart + 0.5 * (det1 - exppart1);
			
			if(j==0){
				loglikeli = lik1;
			}
			else{
				for(int i=0; i<md; i++){
					window[i]=yL_yLM[j+i];
					inp[i]=inputs[j+i];
				}
				
				kernelinv=makeCovarMatrix(inp, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e, &det);
				exppart=(window.transpose() * kernelinv).dot(window);
				lik = md * pipart + 0.5 * (det - exppart);
				
				//conditional likelihood of the last included element
				loglikeli+= (lik1-lik);
			}
		}
	}
	else{
		//full-scale solution
		double det;
		Eigen::MatrixXd kernelinv = makeCovarMatrix(inputs, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e, &det);
		Eigen::VectorXd yL_yLMv (dim);
		for(int i=0; i<dim; i++){
			yL_yLMv[i]=yL_yLM[i];
		}
		double exppart=(yL_yLMv.transpose() * kernelinv).dot(yL_yLMv);
		loglikeli = md * pipart + 0.5 * (det - exppart);
	}

	return -loglikeli; //for minimisation
}

std::vector<std::string> iWQBiasIDARLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	result.push_back("Y_"+varname);
	result.push_back("Ytr_"+varname);
	result.push_back("YB_"+varname);
	result.push_back("YBtr_"+varname);
	result.push_back("YBE_"+varname);
	result.push_back("YBEtr_"+varname);
	result.push_back("I_"+varname);
	
	return result;
}

void iWQBiasIDARLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	int totaldim=mDataTable->numRows();
	
	std::vector<double> Ytrs (totaldim);
	std::vector<double> YBEtrs (totaldim);
	std::vector<double> YBtrs (totaldim);
	
	std::vector<double> Ys (totaldim);
	std::vector<double> YBEs (totaldim);
	std::vector<double> YBs (totaldim);
	
	std::vector<double> Is (totaldim);
	
		
	std::vector<double> past_inputs;
	std::vector<double> future_inputs;
	
	// PART 1: realizations for the past
	// inputs
	inputptr=mDataTable->portForColumn(inputfieldname);
	mDataTable->rewind();
	while(mDataTable->stepRow()!=-1){
		if(mComparisonLink.numeric()){
			past_inputs.push_back(*inputptr);
		}
		else{
			future_inputs.push_back(*inputptr);
		}
	}
	int dim=past_inputs.size();		//size of past only
	
	//residuals
	Eigen::VectorXd yL_yLM (dim);
	mDataTable->rewind();
	int i=0;
	//full length
	while(mDataTable->stepRow()!=-1){
		double model = mComparisonLink.model();
		double modeltr=boxcox_transform(lambda_1, lambda_2, model, NULL);
		if(mComparisonLink.numeric() && i<dim){	//past
			double meas=mComparisonLink.measurement();
			double meastr=boxcox_transform(lambda_1, lambda_2, meas, NULL);
			yL_yLM[i]=meastr-modeltr;
			YBEs[i]=mComparisonLink.measurement();	//this is actually the measurement
			YBEtrs[i]=meastr;
		}
		Ys[i]=model;
		Ytrs[i]=modeltr;
		i++;
	}
	
	//make (E-1 + B-1)-1 in the proper size
	int md;
	int prop_maxkernel = (int)maxkernelsize;
	if(prop_maxkernel % 1 == 0){
		prop_maxkernel++;	//must be odd for inflation
	}
	if(prop_maxkernel < 5){
		prop_maxkernel=5;
	}
	int MAX_KERNEL_SIZE = prop_maxkernel;	
	if(dim<MAX_KERNEL_SIZE){
		md=dim;
	}
	else{
		md=MAX_KERNEL_SIZE;
	}
		
	//full sized realization of B & E for the _past_
	Eigen::MatrixXd SIGMA = inflatedVarBRealization(past_inputs, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e, md);
	Eigen::MatrixXd L=SIGMA.llt().matrixL();
	Eigen::VectorXd indeps (dim);
	for(int i=0; i<dim; i++){
		indeps(i)=invnormdist(0.0, 1.0);
	}
	//multiply SIGMA with SIGMA_E_INV manually (post-multiplication: col-wise)
	for(int r=0; r<dim; r++){
		double invvar = 1.0 / varianceOfE(past_inputs[r], sigma_e2, kappa_e);
		for(int c=0; c<dim; c++){
			SIGMA(c,r) *= invvar;
		}
	}
	
	Eigen::VectorXd mu = SIGMA * yL_yLM; 	//((1.0/sigma_e2) * SIGMA) * yL_yLM;	//was without input dependence in E
	Eigen::VectorXd B (dim);
	Eigen::VectorXd LZ= L * indeps;
	B = mu + LZ;
	
	for(int i=0; i<dim; i++){		//past 
		YBtrs[i]=Ytrs[i] + B[i];
		YBs[i]=boxcox_retransform(lambda_1, lambda_2, YBtrs[i], NULL);
		//Btrs[i]=YBtrs[i]-Ytrs[i];
		if(i==0){
			Is[i]=0.0;
		}
		else{
			double rho = exp(-beta);
			Is[i]=( B[i] - rho * B[i-1]) / sqrt( jumpVarianceOfB(sigma_b2, beta, kappa, pi, past_inputs[i]) );
		}
	}

	//PART 2: make bias & noise process for the _future_
	int future_dim = future_inputs.size();	//size of future
	for(int i=0; i<future_dim; i++){
		double prev_val = YBtrs[dim+i-1]-Ytrs[dim+i-1];
		double jump_var = jumpVarianceOfB(sigma_b2, beta, kappa, pi, future_inputs[i]);
		double newB=makeOUStep(prev_val, jump_var, beta);
		double newE=makeNoiseStep(sigma_e2, future_inputs[i], kappa_e);
		YBEtrs[dim+i]=Ytrs[dim+i]+newB+newE;
		YBEs[dim+i]=boxcox_retransform(lambda_1, lambda_2, YBEtrs[dim+i], NULL);
		YBtrs[dim+i]=Ytrs[dim+i]+newB;
		YBs[dim+i]=boxcox_retransform(lambda_1, lambda_2, YBtrs[dim+i], NULL);
		Is[dim+i]=( newB - exp(-beta)*prev_val) / sqrt( jump_var );
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("Ytr_"+varname)=Ytrs;
	storage->operator[]("YB_"+varname)=YBs;
	storage->operator[]("YBtr_"+varname)=YBtrs;
	storage->operator[]("YBE_"+varname)=YBEs;
	storage->operator[]("YBEtr_"+varname)=YBEtrs;
	storage->operator[]("I_"+varname)=Is;

}

//--------------------------------------------------------------------------------------------------

#pragma mark First-order autoregressive error model with SEP innovations

void iWQARSEPLikelihoodEvaluation::initDefaultParams()
{
	//inits the built-in distribution
	dist.setBeta(0.0);
	dist.setXi(1.0);
	sigma0 = 1.0;
	sigma1 = 0.0;
	mu = 1.0;
	fi = 0.0;
		
	lambda_1 = 1.0;	//defaults to no transform
	lambda_2 = 0.0;
	sumlogy = -DBL_MAX;
}

void iWQARSEPLikelihoodEvaluation::setParams(iWQSettingList list)
{
	std::string varname=mComparisonLink.modelField();
	setParamValueFromMap(&sigma0,"sigma_0",&list, varname);
	setParamValueFromMap(&sigma1,"sigma_1",&list, varname);
	setParamValueFromMap(&beta,"beta",&list, varname);
	setParamValueFromMap(&xi,"xi",&list, varname);
	setParamValueFromMap(&fi,"fi",&list, varname);
	setParamValueFromMap(&mu,"mu",&list, varname);
	setParamValueFromMap(&lambda_1,"lambda_1",&list, varname);	//transformation parameters
	setParamValueFromMap(&lambda_2,"lambda_2",&list, varname);
}

double iWQARSEPLikelihoodEvaluation::evaluate(int startindex, int endindex)
{
	//log likelihood with ARSEP error model
	double loglikeli=0.0;
		
	double prev_bias=0.0;
	
	dist.setBeta(beta);
	dist.setXi(xi);
	
	//get the sum of log values of observations if not done so before
	if(sumlogy==-DBL_MAX){
		sumlogy = 0.0;
		for(int j=startindex; j<endindex; j++){
			mDataTable->setRow(j);
			if(mComparisonLink.numeric()){
				double meas_raw = mComparisonLink.measurement();
				if(meas_raw + lambda_2 > 0.0){
					sumlogy += log(meas_raw + lambda_2);
				}
				else{
					printf("[Warning]: Measurement (%s=%lf at index %d) is not strictly positive after adding lambda_2, so cannot account for lambda_1 in likelihood.\n", measuredFieldName().c_str(), meas_raw, j);
					sumlogy=0.0;
					break;
				}
			}
		}
	}	
		
	loglikeli += (lambda_1 - 1.0) * sumlogy;
	
	for(int j=startindex; j<endindex; j++){
		mDataTable->setRow(j);
		if(mComparisonLink.numeric()){
			double meas=mComparisonLink.measurement();
			double model=mComparisonLink.model();
			double meastr=boxcox_transform(lambda_1, lambda_2, meas, NULL);
			double modeltr=boxcox_transform(lambda_1, lambda_2, model, NULL);
			double act_bias=meastr-modeltr;	
			double innovation = act_bias - fi * prev_bias;
			
			double sigma_t = sigma0 + sigma1 * pow(modeltr > 0.0 ? modeltr : 0.0, mu);
			if(sigma_t <= 0.0){
				sigma_t = 1.0;
			}
			double std_innovation = innovation/sigma_t;
			loglikeli+=dist.logLikeli(std_innovation) - log(sigma_t);
			
			prev_bias=act_bias;
		}
	}
	if(std::isinf(loglikeli) || std::isnan(loglikeli)){
		return DBL_MAX;
	}
	return -loglikeli;	//to make it reversed for minimization
}

std::vector<std::string> iWQARSEPLikelihoodEvaluation::sampleSeriesNames()
{
	std::vector<std::string> result;
	std::string varname=mComparisonLink.modelField();
	result.push_back("Y_"+varname);
	result.push_back("YB_"+varname);
	result.push_back("Ytr_"+varname);
	result.push_back("YBtr_"+varname);
	result.push_back("I_"+varname);
	return result;
}

void iWQARSEPLikelihoodEvaluation::createSampleSeries(std::map<std::string, std::vector<double> > * storage)
{
	if(!storage){
		return;
	}
	
	dist.setBeta(beta);
	dist.setXi(xi);
	
	std::vector<double> Ys;
	std::vector<double> YBs;
	std::vector<double> Ytrs;
	std::vector<double> YBtrs;
	std::vector<double> Is;
	
	mDataTable->rewind();
	double prev_bias=0.0;
	while(mDataTable->stepRow()!=-1){
		double model=mComparisonLink.model();
		double modeltr=boxcox_transform(lambda_1, lambda_2, model, NULL);
		Ys.push_back(model);
		Ytrs.push_back(modeltr);
		if(mComparisonLink.numeric()){
			//past
			double meas=mComparisonLink.measurement();
			double meastr=boxcox_transform(lambda_1, lambda_2, meas, NULL);
			double act_bias=meastr-modeltr;
			YBs.push_back(meas);
			YBtrs.push_back(meastr);
			
			//random increment of the past
			double sigma_t = sigma0 + sigma1 * pow(modeltr > 0.0 ? modeltr : 0.0, mu);
			double condmean = fi * prev_bias;
			double jump = (act_bias - condmean)/(sigma_t!=0.0?sigma_t:1.0);
			Is.push_back(jump);
			
			prev_bias=act_bias;
		}
		else{
			//prediction stage
			double sigma_t = sigma0 + sigma1 * pow(modeltr > 0.0 ? modeltr : 0.0, mu);
			double condmean = fi * prev_bias;
			double val=condmean + sigma_t * dist.generate();
						
			YBs.push_back(boxcox_retransform(lambda_1, lambda_2, modeltr+val, NULL));
			YBtrs.push_back(modeltr+val);
						
			double jump = (val - condmean)/(sigma_t!=0.0?sigma_t:1.0);
			Is.push_back(jump);
			
			prev_bias=val;
		}
	}
	
	std::string varname=mComparisonLink.modelField();
	storage->operator[]("Y_"+varname)=Ys;
	storage->operator[]("YB_"+varname)=YBs;
	storage->operator[]("Ytr_"+varname)=Ytrs;
	storage->operator[]("YBtr_"+varname)=YBtrs;
	storage->operator[]("I_"+varname)=Is;
}

//------------------------------------------------------------------
