/*
 *  spec_model.h
 *  iwq_quality_oxygen
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/DISSOLVED OXYGEN
 *
 */

#ifndef spec_model_h
#define spec_model_h

#include "model.h"

class iWQModel;

//------------------------------------------------------------------------------------------

class IWQ_MODEL_NAME : public iWQModel
{
private:
	// NOTE: declare here every variable, input, flux and parameter as double
	//boundary fluxes
	double F_DO;
	double C_DO;
	double C_DO_sat;
			
	//parameters
	double C_wwtp;
	double C_raw_sewage;
	double C_parasitic;
	double elevation;
	double T_water_min;
			
	//inputs
	double Q_overflow_sewage;
	double Q_overflow_storm;
	double Q_overflow_parasitic;
	double Q_wwtp;
	double T_air;
	double Q_total;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return true; }	
};

//------------------------------------------------------------------------------------------

#endif

