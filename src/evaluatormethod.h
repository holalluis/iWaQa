/*
 *  evaluatormethod.h
 *  Various likelihood calculators
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */
 
#include <map>
#include <string>
#include <vector>

#ifndef evaluatormethod_h
#define evaluatormethod_h

#include "mathutils.h"
#include "complink.h"

class iWQDataTable;
//class iWQCOmparisonLink;
class iWQParameterManager;

typedef std::vector<iWQComparisonLink> iWQComparisonLinkSet;
typedef std::multimap<std::string,std::string> iWQSettingList;	//nested parameter structure

//abstract base class for evaluation methods
class iWQEvaluatorMethod
{
protected:
	iWQDataTable * mDataTable;
	iWQComparisonLink mComparisonLink;			//now all comparison links have their own evaluator method
	iWQParameterManager * mCommonParameters;	//needed for prior likelihood and meta-parameters for error models
	bool setParamValueFromMap(double * dest, std::string key, iWQSettingList * list, std::string flag="");
	bool setParamValueFromMap(std::string * dest, std::string key, iWQSettingList * list, std::string flag="");
		
	//dynamic evaluation parameters
	std::map<std::string, double *> mDynamicParams;	//storage for dynamic evaluation parameters
	bool hasDynamicParamValue(std::string key, double * dest=NULL, std::string flag="");		//responds with true if there is a properly named parameter
	bool setParamValueDynamically(double * dest, std::string key, std::string flag="");	
	
public:
	iWQEvaluatorMethod();
	virtual ~iWQEvaluatorMethod(){ }
	bool wantsParams(){ return (wantsFileParams() || wantsMapParams()); }		
	void setDataTable(iWQDataTable * aTable){ mDataTable=aTable; }
	virtual void setComparisonLink(iWQComparisonLink lnk){ mComparisonLink=lnk; }
	void setParameterStorage(iWQParameterManager * pm){ mCommonParameters=pm; }
	//int dimensions(){ return mComparisonLinks.size(); }
	void updateDynamicParams();	//will be called by the evaluator from outside on parameter update
	double evaluate(); 		//full evaluation by default
	std::string modelFieldName(){ return mComparisonLink.modelField(); }
	std::string measuredFieldName(){ return mComparisonLink.measuredField(); }
	void setLinkPredictiveMode(bool pred){ mComparisonLink.setPredictiveMode(pred); }
		
	//optional implementation
	virtual void setParams(iWQSettingList list){ } 
	virtual void setParams(std::string filename){ }
	virtual bool wantsFileParams(){ return false; }
	virtual bool wantsMapParams(){ return false; }
	virtual void initDefaultParams(){ }
	
	virtual bool isLogScale(){ return false; }	//informs the mcmc sampler about the necessary scaling
	virtual bool priorsApply(){ return false; }	//informs the evaluator to attach priors or not
	
	virtual std::vector<std::string> sampleSeriesNames(){ return std::vector<std::string> (); }	//informs the mcmc sample about the available timeseries
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage){ }			//prepares series samples from the current evaluation
	
	//obligatory implementation
	virtual double evaluate(int startindex, int endindex)=0;
};

//-----------------------------------------------------------------------------------------------

//class factory function
iWQEvaluatorMethod * createEvalMethod(std::string type);

//-----------------------------------------------------------------------------------------------

/*              WARNING:
   FOR TRUE MULTIVARIATE OBJECTIVES 
   we need to use weighing factors in the 
   comparison links, because these 
   methods treat everybody equally. */ 

//-----------------------------------------------------------------------------------------------

//Nash-Sutcliffe index with transform (actually: 1-NS)
class iWQNSBoxCoxEvaluation : public iWQEvaluatorMethod
{
protected:
	double lambda_1;
	double lambda_2;
public:
	virtual void initDefaultParams();	//defaults to no transform
	virtual bool wantsMapParams(){ return true; }
	virtual void setParams(iWQSettingList list); 
	virtual double evaluate(int startindex, int endindex);
};

//-----------------------------------------------------------------------------------------------

//Normal likelihood with transform
class iWQNormalLikelihoodEvaluation : public iWQEvaluatorMethod
{
protected:
	double lambda_1;
	double lambda_2;
	iWQRandomNormalGenerator dist;	//a normal distribution for calculating the likelihood
	double sigma;
	double LOQ;
	double sumlogy;	//sum of log(y_meas_i) to account for lambda_1 in likelihood
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors, no transform
	virtual bool wantsMapParams(){ return true; }
	virtual void setParams(iWQSettingList list); 
	virtual bool isLogScale(){ return true; }
	virtual double evaluate(int startindex, int endindex);
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
	virtual bool priorsApply(){ return true; }
};

//-----------------------------------------------------------------------------------------------

//Heteroscedastic normal likelihood with transform
class iWQHeteroscedasticNormalLikelihoodEvaluation : public iWQNormalLikelihoodEvaluation
{
protected:
	std::string inputfieldname;
	double * inputptr;
	double k_input;
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors, no transform, no driver
	virtual void setParams(iWQSettingList list); 
	virtual double evaluate(int startindex, int endindex);
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
};

//-----------------------------------------------------------------------------------------------

//Normal likelihood with transform on the quantiles of data series
class iWQQuantileNormalLikelihoodEvaluation : public iWQNormalLikelihoodEvaluation
{
private:
	std::vector<double> probs;
	void populateProbs();
protected:
	double quantspacing;						//difference between quantiles (symmetric from 50%)
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors, no transform
	virtual void setParams(iWQSettingList list); 
	virtual double evaluate(int startindex, int endindex);
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
};

//-----------------------------------------------------------------------------------------------

//Normal likelihood with transform on the quantiles of data series
class iWQQuantileLikelihoodEvaluation : public iWQNormalLikelihoodEvaluation
{
private:
	std::vector<double> probs;
	void populateProbs();
protected:
	double quantspacing;						//difference between quantiles (symmetric from 50%)
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors, no transform
	virtual void setParams(iWQSettingList list); 
	virtual double evaluate(int startindex, int endindex);
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
};

//-----------------------------------------------------------------------------------------------

class iWQIDARLikelihoodEvaluation : public iWQEvaluatorMethod
{
//Input-dependent Autoregressive Error model
protected: 
	iWQRandomNormalGenerator dist;				//N(0,1)
	double beta;			
	double sigma_b2;		
	double kappa;			
	double lambda_1;		
	double lambda_2;
	std::string inputfieldname;					//name of the driver
	double * inputptr;
	double sumlogy;	//sum of log(y_meas_i) to account for lambda_1 in likelihood
	
	double transform(double value, bool * error=NULL);
	double retransform(double value, bool * error=NULL);
	
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors
	virtual bool wantsMapParams(){ return true; }
	virtual void setParams(iWQSettingList list); 
	virtual bool isLogScale(){ return true; }
	virtual double evaluate(int startindex, int endindex);	
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
	virtual bool priorsApply(){ return true; }
};

//-----------------------------------------------------------------------------------------------

class iWQBiasIDARLikelihoodEvaluation : public iWQEvaluatorMethod
{
//Input-dependent model bias and indepenedent measurement error
private:
	double maxkernelsize;
protected:
	iWQRandomNormalGenerator dist;	//N(0,1)
	double sigma_b2;				//bias base variance
	double beta;					//log bias correlation
	double kappa;					//bias input dependent variance factor
	double sigma_e2;				//measurement noise variance
	std::string inputfieldname;
	double * inputptr;
	double pi;
	double kappa_e;
	
	//transformation parameters
	double lambda_1;
	double lambda_2;
	double sumlogy;	//sum of log(y_meas_i) to account for lambda_1 in likelihood
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors
	virtual bool wantsMapParams(){ return true; }
	virtual void setParams(iWQSettingList list); 
	virtual bool isLogScale(){ return true; }
	virtual double evaluate(int startindex, int endindex);	
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
	virtual bool priorsApply(){ return true; }
};

//-----------------------------------------------------------------------------------------------

class iWQARSEPLikelihoodEvaluation : public iWQEvaluatorMethod
{
//AR1 process with SEP innovations by Schoups and Vrugt
protected: 
	iWQRandomSEPGenerator dist;	
	double beta;
	double xi;
	double fi;
	double sigma0;
	double sigma1;
	double mu;
	double lambda_1;		
	double lambda_2;
	double sumlogy;	//sum of log(y_meas_i) to account for lambda_1 in likelihood
	
public:
	virtual void initDefaultParams();	//defaults to standard normally distributed errors
	virtual bool wantsMapParams(){ return true; }
	virtual void setParams(iWQSettingList list); 
	virtual bool isLogScale(){ return true; }
	virtual double evaluate(int startindex, int endindex);	
	virtual std::vector<std::string> sampleSeriesNames();
	virtual void createSampleSeries(std::map<std::string, std::vector<double> > * storage);
	virtual bool priorsApply(){ return true; }
};

//-----------------------------------------------------------------------------------------------

#endif
