/*
 *  datatable.cpp
 *  Data table class (like an R data frame) to store and manage data
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/DATA
 *
 */
 
#include "datatable.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#pragma mark string tokenizer

//C++ tokenizer from http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
	tokens.clear();
	
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

//-------------------------------------------------------------------------------------------------

#pragma mark Cross-platform getline

//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\r':
            c = sb->sgetc();
            if(c == '\n')
                sb->sbumpc();
            return is;
        case '\n':
        case EOF:
            return is;
        default:
            t += (char)c;
        }
    }
}

//##################################################################################################

#pragma mark DataTable

iWQDataTable::iWQDataTable()
{
	emptyInit();
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::emptyInit()
{
	mNumRows=0;
	mNumCols=0;
	mActRow=-1;
	mTIndex=-1;
}

//-------------------------------------------------------------------------------------------------

iWQDataTable::~iWQDataTable()
{
	clear();
}

//-------------------------------------------------------------------------------------------------

iWQDataTable::iWQDataTable(iWQDataTable * atable)
{
	emptyInit();
	initFromTable(atable);
}

//-------------------------------------------------------------------------------------------------

iWQDataTable::iWQDataTable(std::string filename)
{
	emptyInit();
	initFromFile(filename);
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::initFromTable(iWQDataTable * atable)
{
	clear();
	//import headers and data
	mColIndexes=atable->mColIndexes;
	mDataStorage=atable->mDataStorage;
	for(int i=0; i<atable->mDataPort.size(); i++){
		double * aval=new double;
		mDataPort.push_back(aval);
	}
	//other attributes
	mNumRows=atable->mNumRows;
	mNumCols=atable->mNumCols;
	mActRow=-1;
	mTIndex=atable->mTIndex;
}

//-------------------------------------------------------------------------------------------------

bool iWQDataTable::isPortValid(double * port)
{
	for(int i=0; i<mNumCols; i++){
		if(mDataPort[i]==port){
			return true;
		}
	}
	return false;
}	

//-------------------------------------------------------------------------------------------------

void iWQDataTable::clear()
{
	for(int i=0; i<mDataPort.size(); i++){
		if(mDataPort[i]){
			delete mDataPort[i];
		}
	}
	for(int i=0; i<mDataStorage.size(); i++){
		mDataStorage[i].clear();
	}
	mNumRows=0;
	mNumCols=0;
	mActRow=-1;
	mTIndex=-1;
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::clearColumn(std::string colname)
{
	int colindex=getColIndex(colname);
	if(colindex!=-1){
		mDataStorage[colindex].assign(mNumRows,0.0);
	}
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::deleteColumn(std::string colname)
{
	int colindex=getColIndex(colname);
	if(colindex!=-1 && colindex!=mTIndex){
		//remove data
		mDataStorage.erase(mDataStorage.begin()+colindex);
		//remove port
		if(colindex<mDataPort.size()){
			if(mDataPort[colindex]){
				delete mDataPort[colindex];
			}
			mDataPort.erase(mDataPort.begin()+colindex);
		}
		else{
			printf("[Error]: Index for column \"%s\" (%d) is out of bounds (%zd).\n",colname.c_str(),colindex,mDataPort.size());
		}
		//remove name
		mColIndexes.erase(colname);
		//reduce col count
		mNumCols--;
	}
	if(colindex!=-1 && colindex==mTIndex){
		printf("[Error]: Cannot delete time column \"%s\".\n",colname.c_str());
	}
}
	
//-------------------------------------------------------------------------------------------------

void iWQDataTable::initFromFile(std::string filename)
{	
	std::ifstream f;
	f.open(filename.c_str(), std::ios::in);
	if(!f.is_open()){
		printf("[Error]: failed to open data file \"%s\".\n",filename.c_str());
		return;
	}
	
	//get rid of previous content
	clear();
	
	if(f.eof()){
		//empty file
		return;
	}
	//read header
	std::string headerline;
	std::getline(f,headerline);
	//safeGetline(f, headerline);
	//explode columns headers
	std::stringstream strstr (headerline);
	std::istream_iterator<std::string> it (strstr);
	std::istream_iterator<std::string> end;
	std::vector<std::string> tokens (it, end);
	
	//for each column make storage, register them in the index and open a port
	for(int i=0; i<tokens.size(); i++){
		addColumn(tokens[i]);
	}
	
	//now read in the content
	std::string dataline;
	std::vector<double> values;
	std::vector<std::string> svalues;
	int linenumber=0;
	while(!f.eof()){
		//read the next data row
		std::getline(f,dataline);
		if(dataline.size()){
			values.clear();
			linenumber++;
			Tokenize(dataline, svalues, " \t");
			if(svalues.size()==mDataStorage.size()){
				//valid number of items in row	
				for(int i=0; i<svalues.size(); i++){
					char * sptr;
					double value=strtod(svalues[i].c_str(), &sptr);
					if(sptr!=svalues[i].c_str()){
						values.push_back(value);
					}
					else{
						//if not valid number, don't complain, just load NaN
						values.push_back(iWQNaN);
					}
				}
				//load them
				mNumRows++;
				for(int i=0; i<mNumCols; i++){
					mDataStorage[i].push_back(values[i]);
				}
			}
			else{
				printf("[Warning]: Row %d contains %zd columns instead of %zd - skipped.\n",linenumber,svalues.size(),mDataStorage.size());
			}
		}
	}	
	f.close();
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::reloadFromFile(std::string filename)
{
	//same as initFromFile but with keeping the structure
	std::ifstream f;
	f.open(filename.c_str(), std::ios::in);
	if(!f.is_open()){
		printf("[Error]: failed to open data file \"%s\".\n",filename.c_str());
		return;
	}
	
	//get rid of previous content
	for(int i=0; i<mDataStorage.size(); i++){
		mDataStorage[i].clear();
	}
	mNumRows=0;
	mActRow=-1;
	
	if(f.eof()){
		//empty file
		return;
	}
	//read header
	std::string headerline;
	std::getline(f,headerline);
	//explode columns headers
	std::stringstream strstr (headerline);
	std::istream_iterator<std::string> it (strstr);
	std::istream_iterator<std::string> end;
	std::vector<std::string> tokens (it, end);
	
	std::vector<int> fieldIndexInOriginal;
	
	
	//for each column make storage, register them in the index and open a port
	for(int i=0; i<tokens.size(); i++){
		fieldIndexInOriginal.push_back(getColIndex(tokens[i]));
	}
	
	//now read in the content
	std::string dataline;
	std::vector<double> values;
	std::vector<std::string> svalues;
	int linenumber=0;
	while(!f.eof()){
		//read the next data row
		std::getline(f,dataline);
		if(dataline.size()){
			values.clear();
			linenumber++;
			Tokenize(dataline, svalues, " \t");
			//int skipped=0;
			if(svalues.size()==mDataStorage.size()){
				//valid number of items in row	
				for(int i=0; i<svalues.size(); i++){
					char * sptr;
					double value=strtod(svalues[i].c_str(), &sptr);
					if(sptr!=svalues[i].c_str()){
						values.push_back(value);
					}
					else{
						//if not valid number, don't complain, just load NaN
						values.push_back(iWQNaN);
					}
				}
				//load them
				mNumRows++;
				for(int i=0; i<mNumCols; i++){
					mDataStorage[fieldIndexInOriginal[i]].push_back(values[i]);
				}
			}
			else{
				printf("[Warning]: Row %d contains %zd columns instead of %zd - skipped.\n",linenumber,svalues.size(),mDataStorage.size());
			}
		}
	}	
	f.close();
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::addColumn(std::string colname, bool warn)
{
	//check if this name already exists
	if(getColIndex(colname)!=-1){
        if(warn){
            printf("[Warning]: Column %s already exists.\n", colname.c_str());
        }
		return;
	}
    std::vector<double> data;
	double * port=new double;
	int index=mNumCols;
	data.assign(mNumRows,0.0);
	mDataStorage.push_back(data);
	mDataPort.push_back(port);
	mColIndexes[colname]=index;
	mNumCols++;
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::addColumns(std::vector<std::string> colnames, bool warn)
{
    for(int i=0; i<colnames.size(); i++){
        addColumn(colnames[i], warn);
    }
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::copyColumn(std::string origin, std::string destination, bool warn)
{
	int srccol=getColIndex(origin);
	if(srccol==-1){
		printf("[Warning]: Column %s does not exist.\n", origin.c_str());
		return;
	}
	
	int destcol=getColIndex(destination);
	if(destcol!=-1 && warn){
		printf("[Warning]: Column %s will be overwritten.\n", destination.c_str());
		return;
	}
	else{
		addColumn(destination, warn);
		destcol=getColIndex(destination);
	}
	if(destcol!=-1){
		//copy data
		mDataStorage[destcol].assign(mDataStorage[srccol].begin(),mDataStorage[srccol].end());
	}
	else{
		printf("[Warning]: Failed to create column %s.\n", destination.c_str());
	}
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::addRows(int count)
{
	if(count<=0){
		return;
	}
	std::vector<double> empty (count,0.0);
	for(int i=0; i<mNumCols; i++){
		mDataStorage[i].insert(mDataStorage[i].end(), empty.begin(), empty.end());
	}
	mNumRows+=count;
}

//-------------------------------------------------------------------------------------------------

int iWQDataTable::numRows()
{
	return mNumRows;
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::setRow(int index)
{
	//is it the previous index?
	if(index==mActRow){
		return;
	}
	commit();
	
	if(index>=0 && index<mNumRows){
		//read values from the new location    
		for(int i=0; i<mNumCols; i++){
			*(mDataPort[i])=mDataStorage[i][index];
		}
			
		//set the index
		mActRow=index;
	}
	else{
		//override with zeros
		for(int i=0; i<mNumCols; i++){
			*mDataPort[i]=0.0;
		}
		mActRow=-1;
	}
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::commit()
{
	//copy (modified) contents back to the storage
	if(mActRow!=-1){
		for(int i=0; i<mNumCols; i++){
			mDataStorage[i][mActRow]=*(mDataPort[i]);
		}
	}
}

//-------------------------------------------------------------------------------------------------

int iWQDataTable::stepRow()
{
	if(mActRow>=0 && mActRow<mNumRows-1){
		setRow(mActRow+1);
	}
	else{
		setRow(-1);
	}
	return mActRow;
}

//-------------------------------------------------------------------------------------------------

int iWQDataTable::getColIndex(std::string colname)
{
	int index=-1;
	if(mColIndexes.find(colname)!=mColIndexes.end()){
		index=mColIndexes[colname];
	}
	return index;
}

//-------------------------------------------------------------------------------------------------

std::string iWQDataTable::nameForColumn(int colindex)
{
    //find the right index
    std::map<std::string, int>::const_iterator it;
    for(it=mColIndexes.begin(); it!=mColIndexes.end(); ++it){
        if(it->second==colindex){
            return it->first;
        }
    }
    return std::string("");
}

//-------------------------------------------------------------------------------------------------

int iWQDataTable::numCols()
{
    return mNumCols;
}

//-------------------------------------------------------------------------------------------------

double * iWQDataTable::portForColumn(std::string colname)
{
	int index=getColIndex(colname);
	return (index>=0 && index<mNumCols)?mDataPort[index]:NULL;
}

//-------------------------------------------------------------------------------------------------

std::string iWQDataTable::columnForPort(double * port)
{
	for(int i=0; i<mDataPort.size(); i++){
		if(mDataPort[i]==port){
			return colNameForIndex(i);
		}
	}
	return "";
}

//-------------------------------------------------------------------------------------------------

double iWQDataTable::valueForColumn(std::string colname)
{
	double * port=portForColumn(colname);
	return (port!=NULL)?*port:0.0;
}

//-------------------------------------------------------------------------------------------------

double iWQDataTable::valueForColumn(std::string colname, int rowindex)
{
	//get the column index
	int index=getColIndex(colname); 
	if(index>=0 && index<mNumCols && rowindex>=0 && rowindex<mNumRows){
		return mDataStorage[index][rowindex];
	}
	else{
		return 0.0;
	}
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::setValueForColumn(double value, std::string colname)
{
	setValueForColumn(value, colname, mActRow);
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::setValueForColumn(double value, std::string colname, int rowindex)
{
	int index=getColIndex(colname);
	if(index>=0 && index<mNumCols && rowindex>=0 && rowindex<mNumRows){
		mDataStorage[index][rowindex]=value;
		if(rowindex==mActRow){
			//refresh THAT value in DataPort
			*mDataPort[index]=value;
		}
	}
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::addToValueForColumn(double value, std::string colname, int rowindex)
{
    int index=getColIndex(colname);
    if(index>=0 && index<mNumCols && rowindex>=0 && rowindex<mNumRows){
        mDataStorage[index][rowindex]+=value;
        if(rowindex==mActRow){
            //refresh THAT value in DataPort
            *mDataPort[index]=mDataStorage[index][rowindex];
        }
    }
}

//-------------------------------------------------------------------------------------------------

std::string iWQDataTable::colNameForIndex(int index)
{
	std::map<std::string, int>::iterator it;
	for(it=mColIndexes.begin(); it!=mColIndexes.end(); ++it){
		if(it->second==index){
			return it->first;
		}
	}
	return "";
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::writeToFile(std::string filename, std::vector<int> colindicestoprint)
{
	FILE * f=fopen(filename.c_str(),"w");
	if(!f){
		printf("[Error]: Failed to write data table to \"%s\".\n",filename.c_str());
        return;
	}
    writeToFile(f,colindicestoprint);
    fclose(f);
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::writeToFile(FILE * f, std::vector<int> colindicestoprint)
{
    if(!f){
        return;
    }
    //actualise data
    commit();
    //print header
    for(int i=0; i<colindicestoprint.size(); i++){
        std::string act_name=colNameForIndex(colindicestoprint[i]);
        if(i>0){
            fprintf(f,"\t");
        }
        fprintf(f,"%s",act_name.c_str());
    }
    fprintf(f,"\n");
    
    //print data
    for(int j=0; j<mNumRows; j++){
        for(int i=0; i<colindicestoprint.size(); i++){
            int act_index=colindicestoprint[i];
            if(act_index>=0 && act_index<mNumCols){
                if(i>0){
                    fprintf(f,"\t");
                }
                double val=mDataStorage[act_index][j];
                if(isnan(val)){
                    fprintf(f,"NA");
                }
                else{
                    fprintf(f,"%0.9lg",val);
                }
            }
        }
        fprintf(f,"\n");
    }
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::saveUNCSIMFormatToFile(std::string filename, std::vector<int> colindicestoprint)
{
	if(mTIndex<0 || mTIndex>=mNumCols){
		printf("[Error]: Cannot export into UNCSIM format without a previously specified time field.\n");
		return;
	}
	FILE * f=fopen(filename.c_str(),"w");
	if(!f){
		printf("[Error]: Failed to write data table to \"%s\".\n",filename.c_str());
		return;
	}
	//actualise data
	commit();
	//print data
	for(int i=0; i<colindicestoprint.size(); i++){
		int index=colindicestoprint[i];
		if(index!=mTIndex){
			std::string act_name=colNameForIndex(index);
			
			std::vector<std::string> data=UNCSIMdata(act_name);			
			for(int j=0; j<data.size(); j++){
				fprintf(f, "%s\n", data[j].c_str());
			}
		}
		//next var comes
	}
	fclose(f);
}

//-------------------------------------------------------------------------------------------------

std::vector<std::string> iWQDataTable::UNCSIMdata(std::string colname, std::string alias)
{
	// Unified UNCSIM coding routine (index coding, no problems with timestamp precision)
	std::vector<std::string> result;
	int index=getColIndex(colname);
	const char * varname=(alias.size()?alias.c_str():colname.c_str());
	char buf [strlen(varname)+20];
	if(index!=-1 && index!=mTIndex){
		for(int j=0; j<mNumRows; j++){
			double val=mDataStorage[index][j];
			if(isnan(val)){
				//silently omit NaNs from the dataset
				//sprintf(buf, "%s_%d\t%0.9lg",varname,j,val);
			}
			else{
				sprintf(buf, "%s_%d\t%0.9lg",varname,j,val);
			}
			result.push_back(buf);
		}
	}
	return result;
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::writeToFile(std::string filename)
{
	writeToFile(filename, getAllIndexes());
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::saveUNCSIMFormatToFile(std::string filename)
{
	saveUNCSIMFormatToFile(filename, getAllIndexes());
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::writeToFile(std::string filename, std::vector<std::string> columnstoprint)
{
	std::vector<int> indexestoprint=getIndexesForColNames(columnstoprint);
	writeToFile(filename, indexestoprint);
}	

//-------------------------------------------------------------------------------------------------

void iWQDataTable::saveUNCSIMFormatToFile(std::string filename, std::vector<std::string> columnstoprint)
{
	std::vector<int> indexestoprint=getIndexesForColNames(columnstoprint);
	saveUNCSIMFormatToFile(filename, indexestoprint);
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::setTField(std::string colname)
{
	mTIndex=getColIndex(colname);
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::setTField(int colindex)
{
	if(colindex>=0 && colindex<mNumCols){
		mTIndex=colindex;
	}
}

//-------------------------------------------------------------------------------------------------

double * iWQDataTable::timePort()
{
	return (mTIndex!=-1)?mDataPort[mTIndex]:NULL;
}

//-------------------------------------------------------------------------------------------------

std::string iWQDataTable::timeColumn()
{
	return colNameForIndex(mTIndex);
}

//-------------------------------------------------------------------------------------------------

std::vector<int> iWQDataTable::getAllIndexes()
{
	std::vector<int> allindices;
	for(int i=0; i<mNumCols; i++){
		allindices.push_back(i);
	}
	return allindices;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> iWQDataTable::getIndexesForColNames(std::vector<std::string> colnames)
{
	std::vector<int> indexestoprint;
	if(mTIndex!=-1){
		indexestoprint.push_back(mTIndex);
	}
	for(int i=0; i<colnames.size(); i++){
		int act_index=getColIndex(colnames[i]);
		if(act_index!=-1 && act_index!=mTIndex){
			indexestoprint.push_back(act_index);
		}
	}
	return indexestoprint;
}

//-------------------------------------------------------------------------------------------------

const std::vector<double> * iWQDataTable::vectorForColumn(std::string colname)
{
	int index=getColIndex(colname);
	if(index!=-1){
		return &(mDataStorage[index]);
	}
	return NULL;
}

//-------------------------------------------------------------------------------------------------

std::vector<std::string> iWQDataTable::columnNames()
{
	std::vector<std::string> result;
	for(int i=0; i<mNumCols; i++){
		result.push_back(colNameForIndex(i));
	}
	return result;
}

//-------------------------------------------------------------------------------------------------

bool iWQDataTable::isRowComplete()
{
	//shows if the current row has NaN(s)
	if(mActRow<0){
		return false;
	}
	bool complete=true;
	for(int i=0; i<mDataStorage.size(); i++){
		if(isnan(mDataStorage[i][mActRow])){
			complete=false;
			break;
		}
	}
	return complete;
}

//-------------------------------------------------------------------------------------------------

#pragma mark Fast searching

void iWQDataTable::createIndexForColumn(std::string colname)
{
    //check if colname is valid
    if(getColIndex(colname)==-1){
        printf("[Error]: Cannot make an index for column %s because it doesn\'t exist.\n",colname.c_str());
        return;
    }
    //check if we already have an index
    if(mValueIndices.find(colname)!=mValueIndices.end()){
        return;
    }
    //make index
    const std::vector<double> * x = vectorForColumn(colname);
    std::vector<int> y (x->size());
    std::size_t n (0);
    std::generate(std::begin(y), std::end(y), [&]{ return n++; });
    std::sort(  std::begin(y),
                std::end(y),
                [&](int i1, int i2) { return (x->at(i1) < x->at(i2)); }
              );
    mValueIndices[colname] = y;
    
    //check if values in the column are unique (need to consider that when searching)
    mValueIndexUnique[colname]=true;
    std::vector<double> z (*x);
    std::sort(std::begin(z), std::end(z));
    for(int i=1; i<z.size(); i++){
        if(z[i]==z[i-1]){
            mValueIndexUnique[colname]=false;
            break;
        }
    }
}

//-------------------------------------------------------------------------------------------------

int iWQDataTable::indexOfKeyValueInColumn(double value, std::string colname)
{
    //find occurrence of a unique value in a specific column
    const std::vector<double> * vals = vectorForColumn(colname);
    if(!vals){
        return -1;
    }
    size_t nrows=vals->size();
    if(nrows==0){
        return -1;
    }
    if(mValueIndices.find(colname)==mValueIndices.end() || mValueIndexUnique[colname]==false){
        // Slow & primitive search algorithm
        for(int i=0; i<nrows; i++){
            if(vals->at(i)==value){
                return i;
            }
        }
        return -1;
    }else{
        //there is an index for this column
        //handle non-existing values
        if(value<vals->at(mValueIndices[colname][0]) || value>vals->at(mValueIndices[colname][nrows-1])){
            return -1;
        }
        if(value==vals->at(mValueIndices[colname][0])){
            return mValueIndices[colname][0];
        }
        std::vector<int>::const_iterator ii = std::lower_bound(mValueIndices[colname].begin(), mValueIndices[colname].end(), value, [&](int i1, double d1){ return ( vals->at(i1) < d1); } );
        long int i = ii - mValueIndices[colname].begin();
        return (int)i;
    }
    return -1;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> iWQDataTable::indexOfValueInColumn(double value, std::string colname)
{
    //find all occurrences of a non-unique value in column
    std::vector<int> result;
    
    //find occurrence of a unique value in a specific column
    const std::vector<double> * vals = vectorForColumn(colname);
    if(!vals){
        return result;
    }
    size_t nrows=vals->size();
    if(nrows==0){
        return result;
    }
    if(mValueIndices.find(colname)==mValueIndices.end()){
        // Slow & primitive search algorithm
        for(int i=0; i<nrows; i++){
            if(vals->at(i)==value){
                result.push_back(i);
            }
        }
        return result;
    }else{
        if(value==vals->at(mValueIndices[colname][0])){
            result.push_back(mValueIndices[colname][0]);
        }
        std::vector<int>::const_iterator iil = std::lower_bound(mValueIndices[colname].begin(), mValueIndices[colname].end(), value, [&](int i1, double d1){ return ( vals->at(i1) < d1); } );
        std::vector<int>::const_iterator iiu = std::upper_bound(mValueIndices[colname].begin(), mValueIndices[colname].end(), value, [&](double d1, int i1){ return ( vals->at(i1) > d1); } );
        
        long int il = iil - mValueIndices[colname].begin();
        long int ul = iiu - mValueIndices[colname].begin();
        
        for(long int i=il; i<ul; i++){
            result.push_back(mValueIndices[colname][i]);
        }
    }
    
    return result;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> iWQDataTable::indicesOfKeyValuesInColumn(std::vector<double> values, std::string colname)
{
    std::vector<int> result;
    size_t nvalues = values.size();
    
    if(mValueIndexUnique[colname]){
        for(int i=0; i<nvalues; i++){
            result.push_back(indexOfKeyValueInColumn(values[i], colname));
        }
    }
    return result;
}

//-------------------------------------------------------------------------------------------------

std::vector<int> iWQDataTable::indicesOfValuesInColumn(std::vector<double> values, std::string colname)
{
    std::vector<int> result;
    size_t nvalues = values.size();
    
    for(int i=0; i<nvalues; i++){
        auto subresult = indexOfValueInColumn(values[i], colname);
        result.insert(result.begin(), subresult.begin(), subresult.end());
    }
    return result;
}

//-------------------------------------------------------------------------------------------------

std::vector<double> iWQDataTable::valuesForIndicesInColumn(std::vector<int> indices, std::string colname)
{
    std::vector<double> result;
    size_t nindices = indices.size();
    for(int i=0; i<nindices; i++){
        result.push_back(valueForColumn(colname, indices[i]));
    }
    return result;
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::createIndexColumnByMatchingColumns(std::string newcolname, std::string colname1, std::string colname2)
{
    //the index of "colname1" where it matches "colname2" will be stored in newcolname
    //one colname1->many colname2 (= colname1 has unique values, colname2 not necessarily)
    //example: newcolname="index_of_neighbour", colname1="id", colname2="neighbour_id"
    addColumn(newcolname);
    
    const std::vector<double> * vals1 = vectorForColumn(colname1);
    const std::vector<double> * vals2 = vectorForColumn(colname2);
    int inew = getColIndex(newcolname);
    std::vector<double> * vnew = &(mDataStorage[inew]);
    if(!vals1 || !vals2){
        return;
    }
    size_t nrows=vals1->size();
    if(nrows==0){
        return;
    }
    
    //now find the index for each value in colname2
    for(int i=0; i<nrows; i++){
        double actval2 = vals2->at(i);
        double d1 = indexOfKeyValueInColumn(actval2, colname1);
        (*vnew)[i]=d1;
    }
}

//-------------------------------------------------------------------------------------------------

void iWQDataTable::print()
{
    printf("\n");
    writeToFile(stdout, getAllIndexes());
    printf("\n");
}


