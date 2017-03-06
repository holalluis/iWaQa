/*
 *  filter.cpp
 *  Time-series filters
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */ 
 
#include <math.h>
#include <stdlib.h>
#include <map>
#include <stdio.h>

#include "filter.h"
#include "datatable.h"

//--------------------------------------------------------------------------------------------------

iWQFilter::iWQFilter()
{
	//provide default functionality: simple copying
	mDataTable = NULL;
	mSrcFieldName = "";
	mDestFieldName = "";
	mFuncName = "copy";
	mAggrFunc = NULL;
	mWindowLength = 1;
	mWindowCenter = 0;
	mSrcCol = NULL;
	mDestPtr = NULL;
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setSrcFieldName(std::string name)
{
	mSrcFieldName = "";
	mSrcCol = NULL;
	if(mDataTable && mDataTable->hasColumnWithName(name)){
		mSrcFieldName = name;
		mSrcCol = mDataTable->vectorForColumn(name);
		return true;
	}
	else{
		if(!mDataTable){
			printf("[Error]: Cannot set filter source field: no data table was specified.\n");
			return false;
		}
		else{
			printf("[Error]: Cannot set filter source field: data table has no field called \"%s\".\n",name.c_str());
			return false;
		}
	}
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setDestFieldName(std::string name)
{
	mDestFieldName = "";
	mDestPtr = NULL;
	if(mDataTable && mDataTable->hasColumnWithName(name)){
		mDestFieldName = name;
		mDestPtr = mDataTable->portForColumn(mDestFieldName);
		return true;
	}
	else{
		if(!mDataTable){
			printf("[Error]: Cannot set filter destination field: no data table was specified.\n");
			return false;
		}
		else{
			printf("[Error]: Cannot set filter destination field: data table has no field called \"%s\".\n",name.c_str());
			return false;
		}
	}
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setDataTable(iWQDataTable * table)
{
	if(!table){
		printf("[Warning]: Setting the data table reference of filter to NULL.\n");
	}
	mDataTable = table;
	return true;
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setFunction(std::string funcname)
{
	mFuncName = "copy";
	mAggrFunc = NULL;
	if(funcname.compare("average")==0){
		mFuncName = funcname;
		mAggrFunc = average;
	}
	if(funcname.compare("variance")==0){
		mFuncName = funcname;
		mAggrFunc = variance;
	}
	if(funcname.compare("min")==0){
		mFuncName = funcname;
		mAggrFunc = min;
	}
	if(funcname.compare("max")==0){
		mFuncName = funcname;
		mAggrFunc = max;
	}
	if(funcname.compare("sumsquares")==0){
		mFuncName = funcname;
		mAggrFunc = sumsquares;
	}
	if(funcname.compare("sum")==0){
		mFuncName = funcname;
		mAggrFunc = sum;
	}
	if(!mAggrFunc){
		printf("[Error]: Cannot set filter type to \"%s\": unknown function. Reverting to default (%s).\n",funcname.c_str(),mFuncName.c_str());
		return false;
	}
	return true;
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setWindowLength(int len)
{
	if(len>=1){
		mWindowLength = len;
		return true;
	}
	else{
		mWindowLength = 1;
		printf("[Error]: Filter window length cannot be smaller than 1 (%d specified).\n",len);
		return false;
	}
}

//--------------------------------------------------------------------------------------------------

bool iWQFilter::setWindowCenter(int center)
{
	if(center>=0){
		mWindowCenter = center;
		return true;
	}
	else{
		mWindowCenter = 0;
		printf("[Error]: Filter window center index cannot be smaller than 0 (%d specified).\n",center);
		return false;
	}
}


//--------------------------------------------------------------------------------------------------

bool iWQFilter::valid()
{
	return (mWindowCenter < mWindowLength && mDataTable && mSrcFieldName.size() && mDestFieldName.size());
}

//--------------------------------------------------------------------------------------------------

void iWQFilter::filter()
{
	//does the actual filtering
	//for each row of destField assemble the vector of sources
	mDataTable->rewind();
	int nrows = mDataTable->numRows();
	for(int r=0; r<nrows; r++){
		mDataTable->setRow(r);	//r is at the center of the window
		int vecstart = r - mWindowCenter;
		int vecend = r + (mWindowLength - mWindowCenter);	//index AFTER the last element
		if(vecstart<0){ vecstart=0; }
		if(vecend>nrows){ vecend=nrows; }
		iWQVector subdata;
		//subdata.assign(mSrcCol->begin() + vecstart, mSrcCol->begin() + vecend);
		for(int s=vecstart; s<vecend; s++){
			subdata.push_back(mSrcCol->at(s));
		}
		double result = mAggrFunc(&subdata);
		*mDestPtr = result;
	}
	mDataTable->commit();
}
