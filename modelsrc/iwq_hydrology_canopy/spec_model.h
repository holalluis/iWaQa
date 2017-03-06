/*
 *  spec_model.h
 *  iwq_hydrology_canopy
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/CANOPY
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
	double intercept_storage;
	
	//boundary fluxes
	double throughfall_mm;
	double et_mm;
	double throughfall_m3s;
	double et_m3s;
	
	//parameters
	double harvest_eff;
	double storage_size;
	double LAI_min;
	double petMult;
	double area;
	
	//inputs
	double rain_mm;
	double pet_mm;

public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
};

//------------------------------------------------------------------------------------------

#endif

