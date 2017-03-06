/*
 *  datatable.h
 *  Data table class (like an R data frame) to store and manage data
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/DATA
 *
 */
 
#include <string> 
#include <vector>
#include <map>
#include <math.h>
#include <limits>

#ifndef datatable_h
#define datatable_h

//--------------------------------------------------------------------------------

//Utilities to have NaN in data tables
#define iWQNaN std::numeric_limits<double>::quiet_NaN()	
 
//--------------------------------------------------------------------------------
 
class iWQDataTable
{
private:
	std::map<std::string, int> mColIndexes;
	std::vector< std::vector<double> > mDataStorage;
	int mNumRows;
	int mNumCols;
	std::vector<double *> mDataPort;
	int mActRow;
	int mTIndex;
	
	int getColIndex(std::string colname);
	std::string colNameForIndex(int index);
	void writeToFile(std::string filename, std::vector<int> colindicestoprint);
    void writeToFile(FILE * f, std::vector<int> colindicestoprint);
	void saveUNCSIMFormatToFile(std::string filename, std::vector<int> colindicestoprint);
		
	std::vector<int> getAllIndexes();
	std::vector<int> getIndexesForColNames(std::vector<std::string> colnames);
	void emptyInit();
    
    std::map<std::string, std::vector<int> > mValueIndices;
    std::map<std::string, bool> mValueIndexUnique;
	
public:
	iWQDataTable();
	iWQDataTable (std::string filename);
	~iWQDataTable();
	iWQDataTable(iWQDataTable * atable);
	void initFromTable(iWQDataTable * atable);
	void initFromFile(std::string filename);
	void reloadFromFile(std::string filename);
	
	void addColumn(std::string colname, bool warn=true);
    void addColumns(std::vector<std::string> colnames, bool warn=true);
	void addRows(int count);
	int numRows();
    int numCols();
    std::string nameForColumn(int colindex);
	void copyColumn(std::string origin, std::string destination, bool warn=true);
	
	void setRow(int index);
	int stepRow();
	int pos(){ return mActRow; }
	void rewind(){ setRow(0); }
	
	void commit();
	
	void clear();
	void clearColumn(std::string colname);
	void deleteColumn(std::string colname);
		
	double * portForColumn(std::string colname);
	std::string columnForPort(double * port);
	double * operator[](std::string colname){ return portForColumn(colname); }
	bool isPortValid(double * port);
	
	void setTField(std::string colname);
	void setTField(int colindex);
	double * timePort();
	std::string timeColumn();
	
	double valueForColumn(std::string colname);
	double valueForColumn(std::string colname, int rowindex);
	void setValueForColumn(double value, std::string colname);
	void setValueForColumn(double value, std::string colname, int rowindex);
    void addToValueForColumn(double value, std::string colname, int rowindex);
	
	void writeToFile(std::string filename);
	void writeToFile(std::string filename, std::vector<std::string> columnstoprint);
	void saveUNCSIMFormatToFile(std::string filename);
	void saveUNCSIMFormatToFile(std::string filename, std::vector<std::string> columnstoprint);
	
	std::vector<std::string> UNCSIMdata(std::string colname, std::string alias="");
	
	const std::vector<double> * vectorForColumn(std::string colname);
	std::vector<std::string> columnNames();
	
	//checking for NaN in data
	bool isRowComplete();					//shows if the current row has NaN(s)
	
	bool hasColumnWithName(std::string colname){ return (portForColumn(colname)!=NULL); }
    
    //seeking for unique key values in a column
    void createIndexForColumn(std::string colname);
    int indexOfKeyValueInColumn(double value, std::string colname);
    std::vector<int> indexOfValueInColumn(double value, std::string colname);
    std::vector<int> indicesOfKeyValuesInColumn(std::vector<double> values, std::string colname);
    std::vector<int> indicesOfValuesInColumn(std::vector<double> values, std::string colname);
    std::vector<double> valuesForIndicesInColumn(std::vector<int> indices, std::string colname);
    
    //printing to stdout
    void print();
    
    //assisting linked hierarchies
    void createIndexColumnByMatchingColumns(std::string newcolname, std::string colname1, std::string colname2);
};

//--------------------------------------------------------------------------------

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");

#endif
