/*
 *  complink.h
 *  Comparison link class 
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */
 
#include <string>
#include <vector>

#ifndef complink_h
#define complink_h

class iWQDataTable;

// Basic link between 2 data fields, a calculated and a measured one 
class iWQComparisonLink
{
private:
	double * modelptr;
	double * measptr;
	std::string modelname;
	std::string measname;
	bool predmode;
public:
	double model();
	double measurement();
	iWQComparisonLink();
	iWQComparisonLink(iWQDataTable * table, std::string modelfield, std::string measfield);
	bool valid(){ return (modelptr && measptr); }
	bool appliesTo(double * port){ return (port==measptr || port==modelptr); }
	double * modelPtr(){ return modelptr; }
	double * measurementPtr(){ return measptr; }
	std::string modelField(){ return modelname; }
	std::string measuredField(){ return measname; }
	bool numeric();					//shows if both numbers are valid
	bool operator==(const iWQComparisonLink & alink) const;
	bool operator!=(const iWQComparisonLink & alink) const;
	void setPredictiveMode(bool p){ predmode=p; }
	bool predictiveMode(){ return predmode; }
};

typedef std::vector<iWQComparisonLink> iWQComparisonLinkSet;

#endif
