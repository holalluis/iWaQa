/*
 *  evaluator.h
 *  General likelihood manager
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */
 
#include <map>
#include <string>
#include <vector>

#ifndef evaluator_h
#define evaluator_h

class iWQSolver;
class iWQParameterManager; 
class iWQInitialValues;
class iWQDataTable;
class iWQComparisonLink;
class iWQEvaluatorMethod;
class iWQRandomGenerator;
class iWQFilter;
class iWQScript;

typedef std::vector<iWQComparisonLink> iWQComparisonLinkSet;
typedef std::map<std::string, double> iWQKeyValues;
typedef std::vector<std::string> iWQStrings;

//-----------------------------------------------------------------------------------------------

// Model evaluator and calibration algorithm
class iWQEvaluator
{
private:
	iWQSolver * mSolver;
	iWQParameterManager * mCommonParameters;
	iWQInitialValues * mInitVals;
	iWQDataTable * mDataTable;
	iWQComparisonLinkSet mComparisonLinks;
	std::vector<iWQEvaluatorMethod *> mEvaluatorMethods;
	std::vector<double> mEvaluatorWeights;
	std::vector<iWQFilter *> mFilters;
	std::vector<iWQScript> mPreScripts;
	std::vector<iWQScript> mPostScripts;
	
	void NelderMead(int n, double start[], double xmin[], double *ynewlo, double reqmin, double step[], int konvge, int kcount, int *icount, int *numres, int *ifault );	
	
	//event-based services: interval indices and state buffer
	int mEvaluateStartRow;
	int mEvaluateEndRow; 
	std::map<std::string, iWQKeyValues> mModelState;
	double * mRainColPtr;	//ptr to the column to adjust
	
public:
	iWQEvaluator();
	~iWQEvaluator();
	
	void setDataTable(iWQDataTable * datatable);
	void setComparisonLinks(iWQComparisonLinkSet links); 
	void setInitialValues(iWQInitialValues * inits);
	void setSolver(iWQSolver * solver);
	void setParameters(iWQParameterManager * parameters);
	void setEvaluatorMethods(std::vector<iWQEvaluatorMethod *> methods);
	void setEvaluatorWeights(std::vector<double> weights);
	void setFilters(std::vector<iWQFilter *> filters);
	void setPreScripts(std::vector<iWQScript> pres);
	void setPostScripts(std::vector<iWQScript> posts);
	
	//accessors
	iWQParameterManager * parameters(){ return mCommonParameters; };
	
	//normal parameter optimization
	void calibrate();
	
	double evaluate(std::vector<double> values);
	double evaluate(double * values, int numpars);
	double evaluate();
	
	bool printWarnings;
	bool returnUnstableSolutions;
	
	//event-based parameter adjustment
	void sequentialCalibrateParameters(std::string eventflagfield);
	void sequentialCalibrateByParameter(std::string eventflagfield);
	//event-based input adjustment
	void sequentialCalibrateInputs(std::string eventflagfield, std::string inputfield, iWQRandomGenerator * inputprior);
	
	//temporary storage for PSO parameters: 
	//TODO: create separate optimizer classes
	bool PSOActive;
	int PSOMaxNumRounds;
	int PSOMaxIdleRounds;
	int PSOSwarmSize;
	
	//temporary storage for nelder-mead simplex parameters
	bool NMSActive;
	int NMSMaxNumRounds;
	double NMSTolerance;
	
	//predictive mode switch
	bool setPredictiveMode(bool mode);
};

#endif
