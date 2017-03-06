/*
 *  spec_model.h
 *  iwq_heat_wastewater
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/WASTEWATER
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

	//variables
	double T_air_smooth;
	double T_ww;
	
	//parameters
	double W_smoothing;
	double T_source_ww;
	double K_ww;
	
	//inputs
	double T_air;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
};

//------------------------------------------------------------------------------------------

#endif

