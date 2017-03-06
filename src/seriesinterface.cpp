/*
 *  seriesinterface.cpp
 *  File interface for 2D time-series
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/DATA
 *
 */ 

#include "datatable.h"
#include <stdlib.h>
#include "seriesinterface.h"
 
iWQSeriesInterface::iWQSeriesInterface(iWQDataTable * tbl)
{
	mDataTable=tbl;
}

//---------------------------------------------------------------------------------------------

iWQSeriesInterface::~iWQSeriesInterface()
{
	removeSeriesLinks();
}

//---------------------------------------------------------------------------------------------

iWQDataTable * iWQSeriesInterface::dataTable()
{
	return mDataTable;
}

//---------------------------------------------------------------------------------------------

void iWQSeriesInterface::removeSeriesLinks()
{
	for(int x=0; x<mColumnMappings.size(); x++){
		if(mColumnMappings[x].fileHandle){
			mColumnMappings[x].fileHandle->close();
			delete mColumnMappings[x].fileHandle;
			mColumnMappings[x].fileHandle=NULL;
		}
	}
	mColumnMappings.clear();
}

//---------------------------------------------------------------------------------------------

void iWQSeriesInterface::addSeriesLink(std::string dataname, std::string filename, bool output)
{
	//check if the data column is valid
	if(!mDataTable){
		printf("[Error]: Cannot create series link to \"%s\" without a valid data table.\n", dataname.c_str());
		return;
	}
	else{
		if(!mDataTable->hasColumnWithName(dataname)){
			printf("[Error]: Series link refers to an invalid column \"%s\".\n",dataname.c_str());
			return;
		}
	}
	
	iWQSeriesMap m;
	m.dataColumnName=dataname;
	m.seriesFileName=filename;
	m.output=output;
		
	//check if this was already here
	bool similar=false;
	int numlinks=mColumnMappings.size();
	for(int x=0; x<numlinks; x++){
		if(mColumnMappings[x].dataColumnName.compare(m.dataColumnName)==0 && mColumnMappings[x].seriesFileName.compare(m.seriesFileName)==0 && mColumnMappings[x].output==m.output){
			similar=true;
			break;
		}
	}
	if(similar){
		printf("[Warning]: series link already defined between data column \"%s\" and file \"%s\". Redefinition ignored.\n", dataname.c_str(), filename.c_str());
		return;
	}
	//not similar
	std::ios_base::openmode mode=output?(std::ios_base::out):(std::ios_base::in);
	std::fstream * fh=new std::fstream(filename.c_str(), mode); 
	m.fileHandle=fh;
	if(!fh->is_open()){
		printf("[Error]: failed to open series file \"%s\" for data column \"%s\".\n",filename.c_str(), dataname.c_str());
		delete fh;
		return;
	}
	mColumnMappings.push_back(m);
}

//---------------------------------------------------------------------------------------------

std::vector<double> iWQSeriesInterface::readALine(std::istream * f)
{
	std::vector<double> result;
	std::string dataline;
	if(f->good() && !f->eof()){
		//read the next data row
		std::getline(*f,dataline);
		std::vector<std::string> tokens;
		if(dataline.size()){
			Tokenize(dataline, tokens, " \t");
			for(int i=0; i<tokens.size(); i++){
				char * sptr;
				double value=strtod(tokens[i].c_str(), &sptr);
				if(sptr!=tokens[i].c_str()){
					result.push_back(value);
				}
				else{
					//if not valid number, don't complain, just load NaN
					result.push_back(iWQNaN);
				}
			}
		}
	}
	return result;
}

//---------------------------------------------------------------------------------------------

void iWQSeriesInterface::writeALine(std::ostream * f, const std::vector<double> * v)
{
	if(!v){
		return;
	}
	int numrecords=v->size();
	if(!f->good() || numrecords==0){
		return;
	}
	
	(*f)<<v->at(0);
	for(int x=1; x<numrecords; x++){
		(*f)<<"\t"<<v->at(x);
	}
	(*f)<<std::endl;
}

//---------------------------------------------------------------------------------------------

void iWQSeriesInterface::refreshInputs()
{
	if(!mDataTable){
		return;
	}
	//get all new columns from series files
	int numlinks=mColumnMappings.size();
	for(int x=0; x<numlinks; x++){
		if(!mColumnMappings[x].output){
			//this is an input stream
			double * ptr=mDataTable->portForColumn(mColumnMappings[x].dataColumnName);
			if(ptr){
				std::vector<double> data=readALine(mColumnMappings[x].fileHandle);
				//now copy data into the dataptr (recycle the data if necessary)
				int datasize=data.size();
				int tablesize=mDataTable->numRows();
				for(int y=0; y<tablesize; y++){
					mDataTable->setRow(y);
					*ptr=data[y%datasize];
				}
				mDataTable->commit();	//after the last update
			}
		}
	}
}

//---------------------------------------------------------------------------------------------

void iWQSeriesInterface::refreshOutputs()
{
	//save specified columns from the data file
	if(!mDataTable){
		return;
	}
	//get all new columns from series files
	int numlinks=mColumnMappings.size();
	for(int x=0; x<numlinks; x++){
		if(mColumnMappings[x].output){
			//this is an output stream
			const std::vector<double> * ptr=mDataTable->vectorForColumn(mColumnMappings[x].dataColumnName);
			if(ptr){
				writeALine(mColumnMappings[x].fileHandle, ptr);
			}
		}
	}
}
