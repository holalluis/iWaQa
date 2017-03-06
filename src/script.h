/*
 *  script.h
 *  External model functionality
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 

#include <string>
#include <vector>

#ifndef iwqscripth
#define iwqscripth

class iWQDataTable;
class iWQParameterManager;

//encapsulation for external scripts or models

class iWQScript
{
	private:
		std::string mCommand;
		std::string mExportTableName;
		std::string mImportTableName;
		std::string mExportParametersName;
		unsigned int mOrder;
		int mReturnStatus;
		iWQDataTable * mDataTable;
		iWQParameterManager * mCommonParameters;
		bool mExportTabDelimitedParamaters;
	public:
		iWQScript();
		
		bool execute();
		
		//property interfaces
		std::string commandString(){ return mCommand; }
		void setCommandString(std::string com){ mCommand=com; }
		std::string exportTableName(){ return mExportTableName; }
		void setExportTableName(std::string etn){ mExportTableName=etn; }
		std::string importTableName(){ return mImportTableName; }
		void setImportTableName(std::string itn){ mImportTableName=itn; }
		std::string exportParametersName(){ return mExportParametersName; }
		void setExportParametersName(std::string epn){ mExportParametersName=epn; }
		void setOrder(unsigned int o){ mOrder=o; }
		unsigned int order(){ return mOrder; }
		int returnStatus(){ return mReturnStatus; }
		iWQDataTable * dataTable(){ return mDataTable; }
		void setDataTable(iWQDataTable * dt){ mDataTable=dt; }
		iWQParameterManager * parameterManager(){ return mCommonParameters; }
		void setParameterManager(iWQParameterManager * pm){ mCommonParameters=pm; }
		
		void setExportTabDelimitedParameters(bool flag){ mExportTabDelimitedParamaters=flag; }
		bool exportTabDelimitedParamaters(){ return mExportTabDelimitedParamaters; }
		
		
		//for sorting
		bool operator<(const iWQScript &rhs) const { return mOrder < rhs.mOrder; }
};

#endif
