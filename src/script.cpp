/*
 *  script.cpp
 *  External model functionality
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 

#include "script.h"

#include "datatable.h"
#include "model.h"

#include <stdlib.h>

//---------------------------------------------------------------------------------

iWQScript::iWQScript()
{
	mCommand = "echo Hello World!";
	mExportTableName = "_data.dat";
	mImportTableName = mExportTableName;
	mExportParametersName = "_pars.dat";
	mOrder = 999;
	mReturnStatus = -9;
	mDataTable=NULL;
	mCommonParameters=NULL;
	mExportTabDelimitedParamaters=false;
}

//---------------------------------------------------------------------------------

bool iWQScript::execute()
{
	//check tools
	if(!mDataTable){
		mReturnStatus=-6;
		return false;
	}
	if(!mCommonParameters){
		mReturnStatus=-5;
		return false;
	}
	if(mCommand.size()==0){
		mReturnStatus=-4;
		return false;
	}
	if(mExportTableName.size()==0){
		mReturnStatus=-3;
		return false;
	}
	if(mImportTableName.size()==0){
		mReturnStatus=-2;
		return false;
	}
	if(mExportParametersName.size()==0){
		mReturnStatus=-1;
		return false;
	}
	//save data table and parameters
	mDataTable->writeToFile(mExportTableName);
	mCommonParameters->saveToFile(mExportParametersName, mExportTabDelimitedParamaters);	//tab delimited export
	
	//run script
	mReturnStatus = system(mCommand.c_str());
	
	//reload data table
	mDataTable->reloadFromFile(mImportTableName);
	
	return (mReturnStatus==0);
}

//---------------------------------------------------------------------------------
