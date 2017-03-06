/*
 *  ascutils.h
 *  ESRI ASCII grid utility functions
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/GRID
 *
 */
 
#include <string>

#ifndef ascutils_h
#define ascutils_h
 
struct asc_header
{
	double xll;
	double yll;
	double cellsize;
	int nrows;
	int ncols;
	double nodata_value;
};

//---------------------------------------------------------------------------------------------------------------

class asc_table
{
public:
	std::string displayname;
	std::string id;
	std::string tabledata;
};

//---------------------------------------------------------------------------------------------------------------

class asc_grid
{
protected:
	asc_header hdr;
	bool integerdata;
	
	void expandstorages(asc_header hdr);
public:
	double ** data;	//pointers to row starts
	double * sdata;	//all values in a single row
	
	std::string id;
	std::string kind;
	std::string displayname;
	
	asc_grid ();
	asc_grid (std::string filename);
	asc_grid (asc_grid * sample, bool takedata=false);
	asc_grid (asc_header hdr1, double ** data1);
	asc_grid (asc_header hdr1, double * sdata1);
	~asc_grid();
	void saveToFile(std::string filename);
	void copyData(double ** data1);
	void copyData(acs_grid * grd){ if(grd && grd->sdata){ copyData(grd->sdata); } };
	
	int sdatasize(){ return hdr.ncols * hdr.nrows * sizeof(double); }
	
	//hdr access methods
	int nrows(){ return hdr.nrows; }
	int ncols(){ return hdr.ncols; }
	double xll(){ return hdr.xll; }
	double yll(){ return hdr.yll; }
	double cellsize(){ return hdr.cellsize; }
	double nodata_value(){ return hdr.nodata_value; }
	
	double min();
	double max();
	
	//additions
	double operator[] (unsigned int index);
	asc_grid& operator= (asc_grid const& g);
	void loadFromFile(std::string filename);
	void setValue(double val);
	void setToNodata(){ setValue(nodata_value()); }
	void fillClassified(std::map<int, double> values, asc_grid * class_map);
	void fillClassified(std::map<int, double *> values, asc_grid * class_map);
	double * ptrToCell(int r, int c);
	double * ptrToCoord(double x, double y);
	void addValuesFrom(asc_grid * grd);
	bool valid(int w){ return (sdata && w>=0 && w<hdr.nrows*hdr.ncols && sdata[w]!=hdr.nodata_value); }
	
	asc_header header(){ return hdr; }
	
	void setNodataValue(double val){ hdr.nodata_value=val; }
	
	//integer flags
	void setIntegerType(bool flag){ integerdata=flag; }
	bool isIntegerType(){ return integerdata; }
};

//---------------------------------------------------------------------------------------------------------------

struct grid_coord
{
	int r;
	int c;
};

//---------------------------------------------------------------------------------------------------------------

template <class T> T ** expandGrid(T ** base, asc_header hdr);
template <class T> void collapseGrid(T ** base, asc_header hdr);
FILE * openasc_w(std::string filename, asc_header hdr);
FILE * openasc_r(std::string filename, asc_header &hdr);

//---------------------------------------------------------------------------------------------------------------

struct point
{
	double x;
	double y;
};

//---------------------------------------------------------------------------------------------------------------

point makepoint(double x, double y);

//---------------------------------------------------------------------------------------------------------------

template <class T>
T ** expandGrid(T ** base, asc_header hdr)
{
	base=new T * [hdr.nrows];
	for(int x=0; x<hdr.nrows; x++){
		base[x]=new T [hdr.ncols];
	}
	return base;
}

//---------------------------------------------------------------------------------------------------------------

template <class T>
void collapseGrid(T ** base, asc_header hdr)
{
	if(!base){
		return;
	}
	for(int x=0; x<hdr.nrows; x++){
		if(base[x]){
			delete [] base[x];
		}
	}
	delete [] base;	
}

#endif

