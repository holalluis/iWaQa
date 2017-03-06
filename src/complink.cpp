/*
 *  complink.cpp
 *  Comparison link class 
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */

#include "complink.h"
#include "datatable.h"
		
#pragma mark ComparisonLink

double iWQComparisonLink::model()
{
	if(modelptr){
		return *modelptr;
	}
	else{
		return 0.0;
	}
}

//--------------------------------------------------------------------------------------------------

double iWQComparisonLink::measurement()
{
	if(measptr){
		return *measptr;
	}
	else{
		return 0.0;
	}
}

//--------------------------------------------------------------------------------------------------

iWQComparisonLink::iWQComparisonLink()
{
	modelptr=NULL;
	measptr=NULL;
	modelname="";
	measname="";
	predmode=false;
}

//--------------------------------------------------------------------------------------------------
	
iWQComparisonLink::iWQComparisonLink(iWQDataTable * table, std::string modelfield, std::string measfield)
{
	modelptr=NULL;
	measptr=NULL;
	
	if(!table){
		return;
	}
	
	modelptr=table->portForColumn(modelfield);
	measptr=table->portForColumn(measfield);
	
	if(modelptr){
		modelname=modelfield;
	}
	if(measptr){
		measname=measfield;
	}
	
	predmode=false;
}


//-------------------------------------------------------------------------------------------------

bool iWQComparisonLink::operator==(const iWQComparisonLink & alink) const
{
	return (modelptr==alink.modelptr && measptr==alink.measptr);
}

//-------------------------------------------------------------------------------------------------

bool iWQComparisonLink::operator!=(const iWQComparisonLink & alink) const
{
	return !(*this == alink);
}

//-------------------------------------------------------------------------------------------------

bool iWQComparisonLink::numeric(){ 
	if(predmode){
		return false;	//pretend that link is not complete
	}
	else{
		return valid() && !(isnan(*modelptr) || isnan(*measptr)); 
	}
}

