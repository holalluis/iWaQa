/*
 *  filter.h
 *  Time-series filters
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */ 
 
#include <string>
#include <vector>

#ifndef filter_h
#define filter_h

#include "mathutils.h"

class iWQDataTable;
class iWQParameterManager;

typedef std::multimap<std::string,std::string> iWQSettingList;	//nested parameter structure

//data filter class
class iWQFilter
{
protected:
	iWQDataTable * mDataTable;
	std::string mSrcFieldName;
	std::string mDestFieldName;
	std::string mFuncName;
	double (*mAggrFunc)(const iWQVector *);	//pointer to an aggregating function
	int mWindowLength;						//total length of the aggregation window
	int mWindowCenter;						//relative position of the output cell
	const iWQVector * mSrcCol;
	double * mDestPtr;
public:
	iWQFilter();
	std::string destFieldName(){ return mDestFieldName; }
	std::string srcFieldName(){ return mSrcFieldName; }
	bool setSrcFieldName(std::string name);
	bool setDestFieldName(std::string name);
	bool setDataTable(iWQDataTable * table);
	iWQDataTable * dataTable(){ return mDataTable; }
	bool setFunction(std::string funcname);
	std::string function(){ return mFuncName; }
	bool setWindowLength(int len);
	int windowLength(){ return mWindowLength; }
	bool setWindowCenter(int len);
	int windowCenter(){ return mWindowCenter; }
	
	bool valid();
	void filter();		//does the actual filtering
};

#endif
