/*
 *  spec_model.h
 *  iwq_quality_ph
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/pH
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
	double pH;
	double C_TIC;
	double F_alk;
	double C_alk;
			
	//parameters
	double pH_wwtp;
	double pH_raw_sewage;
	double pH_rain;
	double pH_natural;
	double C_alk_wwtp;
	double C_alk_raw_sewage;
	double C_alk_natural;
			
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

