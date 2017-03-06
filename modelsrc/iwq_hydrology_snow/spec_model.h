/*
 *  spec_model.h
 *  iwq_hydrology_snow
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/SNOW
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
	double snow;
	
	//boundary fluxes
	double rain;
	
	//parameters
	double rMult;
	double tcrit;
	double tmelt;
	double ksnow;
	double preclambda;
	
	//inputs
	double prec;
	double temp;

public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
};

//------------------------------------------------------------------------------------------

#endif

