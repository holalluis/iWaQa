/*
 *  spec_model.h
 *  iwq_heat_network
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/STREAM NETWORK
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
	
	//fluxes
	double T_water;
	double T_water_flux;
	
	//parameters
	double tauperd;
	double T_offset;
	double shaded;
	double LAI_MAX;
		
	//inputs
	double T_air;
	double T_source;
	double T_eq;
	double T_eq_shade;
	double K_model;
	double K_model_shade;
	double Q;
	double LAI;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
};

//------------------------------------------------------------------------------------------

#endif

