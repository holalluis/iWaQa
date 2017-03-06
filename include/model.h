/*
 *  model.h
 *  generic model functionality
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 
 
#include <vector>
#include <map>
#include <string>

#ifndef model_h
#define model_h

//macro shortcuts for definitions 
#define PAR(x) defineParam( &x , #x )
#define VAR(x) defineVariable( &x , #x , false )
#define INP(x) defineInput( &x, #x )
#define BFX(x) defineVariable( &x , #x , true )
#define CVAR(x) defineConcentrationVariable( &x , #x )		//Concentration variable for transport models

class iWQModel;
class iWQInitialValues;
class iWQRandomGenerator;

#include "lsodaintegrator.h"

//################################################################################

#pragma mark Plugin interface v 1.0

#define IWQ_PLUGIN_INTERFACE_VERSION_MAJOR 1
#define IWQ_PLUGIN_INTERFACE_VERSION_MINOR 0

#define IWQ_MODEL_CREATOR_METHOD iWQModelCreate
#define IWQ_MODEL_DESTROY_METHOD iWQModelDestroy
#define IWQ_MODEL_IDENTIFIER_METHOD iWQModelIdentifier
#define IWQ_PLUGIN_MAJOR_VERSION_METHOD iWQPluginMajorVersion
#define IWQ_PLUGIN_MINOR_VERSION_METHOD iWQPluginMinorVersion

typedef iWQModel * (*iWQModelFactoryMethod)();			//IWQ_MODEL_CREATOR_METHOD
typedef void (*iWQModelDestructorMethod)(iWQModel *);	//IWQ_MODEL_DESTROY_METHOD
typedef std::string (*iWQModelIdentifierMethod)();		//IWQ_MODEL_IDENTIFIER_METHOD
typedef int (*iWQPluginVersionMethod)();				//IWQ_PLUGIN_MAJOR_VERSION_METHOD and IWQ_PLUGIN_MINOR_VERSION_METHOD

//################################################################################

#pragma mark Other classes

// Convenience container types
typedef std::vector<std::string> iWQStrings;
typedef std::map<std::string, double> iWQKeyValues;

class iWQLimits
{
public:
	double min;
	double max;
};

//--------------------------------------------------------------------------------------------

// Protocol for parameter handlers
class iWQParameterHandler
{
public:
	virtual void setValueForParam(double value, std::string key)=0;
	virtual double valueForParam(std::string key)=0;
	virtual bool hasValueForParam(std::string key)=0;
	
	//default bypassing to make the use of flags optional
	virtual void setValueForParam(double value, std::string key, std::string flag){ setValueForParam(value, key); }
	virtual double valueForParam(std::string key, std::string flag){ return valueForParam(key); }
	virtual bool hasValueForParam(std::string key, std::string flag){ return hasValueForParam(key); }	
};

//--------------------------------------------------------------------------------------------

// Shared parameter manager
class iWQParameterManager : public iWQParameterHandler
{
private:
	std::vector<double *> mLocalParams;
	std::map<std::string, double *> mParams;
	std::vector<iWQModel *> mBoundClients;
	std::map<std::string, iWQLimits> mLimits;
	
	std::string makeFlaggedStr(std::string key, std::string flag);
	void decomposeFlaggedStr(std::string flaggedstr, std::string & key, std::string & flag);
	
	std::string makeLiteralFlaggedStr(std::string flaggedstr);
	std::string makeFlaggedStr(std::string literalflaggedstr);
	
	//parameter likelihood distributions
	std::map<std::string, iWQRandomGenerator *> mLinkedDistributions;
	std::vector<iWQRandomGenerator *> mOrderedDistributions;	
	
	int indexOfParam(std::string key);
	
public:
	iWQParameterManager();
	~iWQParameterManager();
	
	//attached parameter handlers
	void bindRequest(iWQModel * client);
	void detachRequest(iWQModel * client);
	
	//define (initialize) parameter
	void defineParam(std::string key);
	void initParam(std::string key, double value);
		
	//define (initialize) per flags
	void defineParam(std::string key, std::string flag);
	void initParam(std::string key, std::string flag, double value);
	
	//clear parameters
	void clearAllParams();
	void clearParam(std::string key);
	void clearParam(std::string key, std::string flag);
	
	//access values
	virtual void setValueForParam(double value, std::string key);
	virtual double valueForParam(std::string key);
	virtual bool hasValueForParam(std::string key);
	
	//limits
	void setLimitsForParam(iWQLimits lim, std::string key, std::string domain);
	void clearLimitsForParam(std::string key);
	bool hasLimitsForParam(std::string key);
	iWQLimits limitsForParam(std::string key);
	
	//access values by flags
	virtual void setValueForParam(double value, std::string key, std::string flag);
	virtual double valueForParam(std::string key, std::string flag);
	virtual bool hasValueForParam(std::string key, std::string flag);
	
	//wildcard access for flagged parameters by grid-type models
	std::map<int, double> valuesForParam(std::string key, std::string flag);
		
	//(un)archiving
	void initFromFile(std::string filename);
	void initFromTabDelimitedFile(std::string filename);
	void saveToFile(std::string filename, bool tabdelimited=false);
	
	//interface for unnamed access of locally stored parameters (for calibration robots)
	int numberOfParams();
	std::vector<double> plainValues();
	void setPlainValues(std::vector<double> new_values);
	void setPlainValues(double * new_values, int numvals);
	std::vector<std::string> namesForPlainValues();
	
	//link distributions to specific parameters
	void linkDistributionToParam(iWQRandomGenerator * dist, std::string key, std::string flag="");
	void detachDistributionFromParam(std::string key, std::string flag="");
	iWQRandomGenerator * distributionForParam(std::string key, std::string flag="");
	double logLikelihood(bool report=false);
	double logLikelihoodOfSet(std::vector<double> values);	
};

//-------------------------------------------------------------------------------------------------------------

typedef double * (*iWQVariableDeltaAccessor)(iWQModel * aModel, double * variable);

double * iWQDelta(iWQModel * aModel, double * variable);
double * iWQDiagnosticDelta(iWQModel * aModel, double * variable);

//-------------------------------------------------------------------------------------------------------------

class iWQModel : public iWQParameterHandler
	{
	private:
		std::map<std::string, double *> mParams;
		std::map<double*, double*> mDerivLocations;
		std::vector<double *> mVarLocations;
		std::vector<double *> mInputLocations;
		std::vector<bool> mShouldTakeDelta;
		iWQStrings mVarNames;
		iWQStrings mInputNames;
		
		iWQParameterManager * mParentParameterManager;
		
		//solver things
		double * mA;
		double ** mB;
		double * mC;
		double * mD;
		double ** mF;
		
		//internal solver interface
		void readVariables(double * from, int length);
		void copyDerivatives(double * dest, int length);
			
		//storage for model flags
		iWQStrings mModelFlags;
		
		//model identifier
		std::string mModelId;
		
		//internal diagnostics
		std::map<std::string,bool> mParamInitState;
		std::map<std::string, bool> mInputConnectedState;
		double * mFooVarDeltaContainer;		//to redirect erroneous delta requests
		bool mVarDeltaErrorIndicator;
		iWQVariableDeltaAccessor dFunction;
		
		friend double * iWQDelta(iWQModel * aModel, double * variable);
		friend double * iWQDiagnosticDelta(iWQModel * aModel, double * variable);
		
		
	protected:	
		//definition methods
        void defineVariable(double * var, const char * varname, bool delta);
        void defineInput(double * inp, const char * inputname);
		void defineParam(double * par, const char * parname);
		
		//connection between a variable and its derivative
		double * d(double & var);
		double * D(double * var);	//method for the abstract channel transport module
		
		//boundary flux connector
		double * F(double & var);
		
		//type identifier
		std::string mTypeId;
		
		//different solvers
		bool solve1StepRungeKuttaFehlberg(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps);	//built in
		bool solve1StepLSODA(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps);				//external
		
		friend class iWQLSODAIntegrator;
		
	public:
		iWQModel(std::string type);
		virtual ~iWQModel();
		
		//parameter management
		virtual void setValueForParam(double value, std::string key);
		virtual double valueForParam(std::string key);
		virtual bool hasValueForParam(std::string key);
		void bind(iWQParameterManager * par);
		void detach();
		bool isBound() const;
		iWQParameterManager * sharedManager();
		virtual void updateParameters();	//virtual because of GRID models
		
		//initial value resolver: public from version 4.1.1
        void setInitialValues(iWQInitialValues * initvals);
	
		//name and dimension interfaces
		iWQStrings outputDataHeaders() const;
		iWQStrings inputDataHeaders() const;
		iWQStrings parameters() const;
		int numVariables() const;
		
		//outlets
		virtual double * rwoutlet(std::string name) const;
		virtual const double * routlet(std::string name) const;
		
		//integration
        virtual bool solve1Step(double xvon, double xbis, iWQInitialValues * yvon, double hmin, double eps);
		virtual bool isStatic(){ return false; }	//return true if the model does not need to be integrated (non-differential type models)
		
		//computation methods
		virtual void modelFunction(double x)=0; 
		
		//model identifier
		void setModelId(std::string newid);
		std::string modelId() const;
		
		//parameter evaluation for user-defined constraints
		virtual bool verifyParameters(){ return true; }
		
		//model flags and types
		iWQStrings modelFlags() const;
		void setModelFlags(iWQStrings flags);
		void setModelFlag(std::string flag);
		std::string modelType() const;
		
		//diagnostics
		bool verify();
		
		//made public for DEBUG
		void setStateVariable(std::string name, double value);
		void resetState();
		iWQStrings variableNames() const;
	};

//-----------------------------------------------------------------------------------------------

// Model class for any channel transport schema (CSTR concept)

class iWQGenericChannelTransport : public iWQModel
{
protected:
	double * mConcentrationVariable;

	//boundary fluxes
	double Fout;
		
	//parameters
	double channel_width;
	double channel_depth;
	double channel_length;
	
	//inputs
	double Fin;
	double Fnew;
	double Qout;

	//utility function to define primary variables
	void defineConcentrationVariable( double * var, const char * varname );
	
public:
	iWQGenericChannelTransport (std::string modelname);
	virtual void modelFunction(double x); 
	virtual ~iWQGenericChannelTransport(){ }
	virtual double R(double x){ return 0; }
	virtual bool verifyParameters(){ return true; }
};

//-----------------------------------------------------------------------------------------------

// Associative container for initial values
class iWQInitialValues
{
private:
	iWQKeyValues mDefaultValues; //for any model
	std::map<std::string, iWQKeyValues> mValues; //key is modelid
	
	iWQParameterManager * mSharedParameterManager;
	
public:
	iWQInitialValues();
	
	bool hasDefaultValueForVariable(std::string varname);
	bool hasValueForVariable(std::string varname, std::string modelid);
	double defaultValueForVariable(std::string varname);
	double valueForVariable(std::string varname, std::string modelid);
	void setDefaultValueForVariable(double value, std::string varname);
	void setValueForVariable(double value, std::string varname, std::string modelid);
	iWQKeyValues variablesForId(std::string modelid);
	iWQKeyValues defaultVariables();
	
	void setParameterManager(iWQParameterManager * par);
	iWQParameterManager * parameterManager();
	
	//to get a dump from all the available values
	iWQKeyValues allValues();
};

//-------------------------------------------------------------------------------------------------------------

#endif
