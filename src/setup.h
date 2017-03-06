/*
 *  setup.h
 *  Highest-level object connecting everything and interpreting layouts
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/SETUP
 *
 */ 

#include <string>
#include <vector>
#include <map>
#include "complink.h"

class iWQModel;
class iWQParameterManager;
class iWQInitialValues;
class iWQLink;
class iWQSolver;
class iWQEvaluator;
class iWQEvaluatorMethod;
class iWQDataTable;
class iWQModelFactory;
class TiXmlHandle;
class TiXmlNode;
class TiXmlElement;
class iWQRandomGenerator;
class iWQSeriesInterface;
class iWQFilter;
class iWQScript;

typedef iWQRandomGenerator iWQDistribution;

// Quality parameter to evaluate run/calibration capability
typedef enum { IWQ_NOT_VALID=0, IWQ_VALID_FOR_RUN=1, IWQ_VALID_FOR_CALIBRATE=2 } iWQModelLayoutValidity;

// A wrapper class to keep a model setup together + import routine from xml file
 
class iWQModelLayout
{
private:
	// Import functions
	void loadModels(TiXmlHandle docHandle);
	void loadParameters(TiXmlHandle docHandle);
	void loadData(TiXmlHandle docHandle);
	void loadConnections(TiXmlHandle docHandle);
	void loadFilters(TiXmlHandle docHandle);
	void loadInitVals(TiXmlHandle docHandle);
	void loadComparisonLinks(TiXmlHandle docHandle); 
	iWQEvaluatorMethod * loadEvaluationMethod(TiXmlElement * compareNode, iWQComparisonLink link, std::string methodName);
	void loadDistributions(TiXmlHandle docHandle);
	void loadScripts(TiXmlHandle docHandle);
	void configureSolver(TiXmlHandle docHandle);
	void configureOptimizer(TiXmlHandle docHandle);
		
	bool checkLayoutVersion(TiXmlHandle docHandle);
		
	void storeNodeInMap(TiXmlNode * pParent, std::multimap<std::string, std::string> * container, std::string prefix, int level);
	
	std::vector<iWQModel *> mModels;
	std::vector<iWQLink> mLinks;
	std::vector<iWQLink> mExportLinks;
	std::vector<std::string> mDataColsToExport;
	std::vector<iWQFilter *> mFilters;
	iWQSolver * mSolver;
	iWQParameterManager * mCommonParameters;
	iWQEvaluator * mEvaluator;
	iWQDataTable * mDataTable;
	iWQInitialValues * mInitVals;
	std::vector<iWQComparisonLink> mComparisonLinks;
	iWQModelFactory * mModelFactory;
	std::vector<iWQEvaluatorMethod *> mEvaluatorMethods;
	std::vector<double> mEvaluatorWeights;
	std::map<std::string, iWQDistribution *> mDistributions;	//named distributions for all purposes
	iWQSeriesInterface * mSeriesInterface;
	
	std::vector<iWQScript> mPreScripts;
	std::vector<iWQScript> mPostScripts;
	
	std::string mFilename;
	void printError(std::string errormessage, TiXmlElement * element, int errorlevel=1);
	
	bool runmodel(int * firsterrorrow=NULL, double * firsterrort=NULL);	//core running routine
	
	void saveBestSolutionSoFar();	//helper for MCMC
	
public:
	//constructor/destructor
	iWQModelLayout(std::string filename);
	~iWQModelLayout();
	
	//validator
	iWQModelLayoutValidity validity();
	
	//accessors for internal components
	std::vector<iWQModel *> models(){ return mModels; }
	std::vector<iWQLink> links(){ return mLinks; }
	std::vector<iWQFilter *> filters(){ return mFilters; }
	iWQSolver * solver(){ return mSolver; }
	iWQParameterManager * parameters(){ return mCommonParameters; }
	iWQEvaluator * evaluator(){ return mEvaluator; }
	iWQDataTable * dataTable(){ return mDataTable; }
	iWQInitialValues * initialValues(){ return mInitVals; }
	std::vector<iWQComparisonLink> comparisonLinks(){ return mComparisonLinks; }
	iWQDistribution * distributionForName(std::string name){ return mDistributions[name]; }
	
	//unified run method
	void run();
	
	//evaluator wrapper methods
	void calibrate();
	double evaluate();
	
	//unified sensitivity analysis routines
	void localSensitivityAnalysis(double rel_deviance, std::string target, std::string filename);
	void regionalSensitivityAnalysis(double rel_deviance, std::string target, std::string filename, int numsimulations);
	
	//file I/O wrappers
	void saveParameters(std::string filename, bool tabdelimited=false);
	void saveResults(std::string filename);
	void saveResultsUNCSIM(std::string filename);
	void saveLayoutGraph(std::string dotfilename);
	void loadParameters(std::string filename, bool tabdelimited=false);
	
	//Markov-chain Monte Carlo sampling (now on the evaluator function)
	void MCMC(int numrounds, int burnin, std::string filename, bool loadpropmatrix=false);			//Plain own adaptive Metropolis
	void MCMC_Haario(int numrounds, int burnin, std::string filename);	//Haario's continually adaptive algorithm from MHAdaptive
	void runOnSample(std::string samplefilename, std::string outputfilename);
	void runStandardSeriesOnSample(std::string samplefilename, int desiredrowcount=1000, bool predictivemode=false, bool binary=true);
	void createBestSeries(std::string parbestfilename);
	
	//UNCSIM configurator
	void furnishUNCSIM(std::string dirname);
	
	//misc. utilities 
	std::string filename(){ return mFilename; }   //document name
	
	//print general model info
	void printModelInfo(std::string name);
	
	//diagnosis
	bool verify(); //check for model defects, uninitialized parameters, unconnected inputs, etc
	
};
