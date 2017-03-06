/*
 *  seriesinterface.h
 *  File interface for 2D time-series
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/DATA
 *
 */ 

#include <string> 
#include <vector>
#include <iostream>
#include <fstream>

#ifndef series_interface_h
#define series_interface_h

class iWQDataTable;

//A container for internal storage
class iWQSeriesMap
{
public:
	std::string dataColumnName;
	std::string seriesFileName;
	std::fstream * fileHandle;
	bool output;
};

class iWQSeriesInterface
{
private:
	std::vector<double> readALine(std::istream * f);
	void writeALine(std::ostream * f, const std::vector<double> * v);
protected:
	iWQDataTable * mDataTable;
	std::vector<iWQSeriesMap> mColumnMappings;
public:
	iWQSeriesInterface(iWQDataTable * tbl);
	~iWQSeriesInterface();
	iWQDataTable * dataTable();
	
	void removeSeriesLinks();
	void addSeriesLink(std::string dataname, std::string filename, bool output=false);
	
	void refreshInputs();
	void refreshOutputs();
};

#endif
