/*
 *  spec_model.h
 *  iwq_hydrology_paved
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/PAVED SURFACES
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
	double storage;
	
	//boundary fluxes
	double et;
	double runoff;
	double infiltration;
	
	//parameters
	double area;
	double s;
	double k_s;
	double petMult;
	double k_infiltr;
    double k_impermeable;
    double s_mult;
		
	//inputs
	double rain;
	double pet;

public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
};

//------------------------------------------------------------------------------------------

#endif

