/*
 *  solver.h
 *  General model solver
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 

#include <vector>
#include <map>
#include <string>

#ifndef solver_h
#define solver_h

#include "model.h"

//-------------------------------------------------------------------------------------------------------------

// Link between models and/or data
class iWQLink
	{
	private:
		iWQModel * srcmod;
		iWQModel * destmod; 
		const double * srcptr;
		double * destptr;
		
		//proportional stuff
		bool fixed_proportion;
		double proportion;
		bool keyed_proportion;
		const double * prop_numerator;
		const double * prop_denominator;
		
		friend class iWQSolver;
		
		void zerodest();
		void linkadd();
		
		void init();	//To zero out all member fields before establishment
		
	public:
		void establish(iWQModel * src, std::string srcname, iWQModel * dest, std::string destname); 	//INPUT FROM ANOTHER MODEL
        void establish(const double * src, iWQModel * dest, std::string destname);                      //INPUT FROM OUTSIDE
        void establish(iWQModel * src, std::string srcname, double * dest);								//JUST FOR OUTPUT
		void establish(const double * src, double * dest);		

		//proportional links
		void setFixedProportion(double prop);
		void setKeyedProportion(std::string key);
		
		//COPY DATA
		iWQModel * dependsOn();
		iWQModel * subject();
		double * destinationPort(){ return destptr; } 	//for graph analysis
        void print();   //for debug
		
		bool operator==(const iWQLink & alink) const;
		bool operator!=(const iWQLink & alink) const;
	};

//-----------------------------------------------------------------------------------------------

// Collection of links	
typedef std::vector<iWQLink> iWQLinkSet;

//-----------------------------------------------------------------------------------------------

class iWQDataTable;

// Solver for a set of linked models
class iWQSolver
	{
	private:
		iWQLinkSet mLinks;
		iWQLinkSet mExportLinks;
		std::vector<iWQModel *> mModels;
		iWQLinkSet mInterLinks;
		bool mTreeError;
		double mHmin;
		double mEps;
		
		std::vector<iWQModel *> mFaultyModels;	//storage for models that did not solve properly

	public:
		iWQSolver(iWQLinkSet inputlinks, iWQLinkSet outputlinks);
		iWQSolver(iWQLinkSet alllinks);
		void setMinStepLength(double value);
		void setAccuracy(double value);
		double minStepLength(){ return mHmin; }
		double accuracy(){ return mEps; } 
        bool solve1Step(double xvon, double xbis, iWQInitialValues * yvon);
		bool valid();
		bool saveInitVals(iWQInitialValues * yfrom);
		std::vector<std::string> exportedDataHeaders(iWQDataTable * datatable);
		std::vector<iWQModel *> modelsThatDidNotSolve();
		
		//not 100% tested but seems to work
		std::map<std::string, iWQKeyValues> modelState();
		void setModelState(std::map<std::string, iWQKeyValues> state);
	};

#endif
