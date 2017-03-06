/*
 *  spec_model.h
 *  iwq_heat_soil
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/SOIL
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
	double T_soil;
	double LAI;
	
	//parameters
	//soil
	double K_soil;
	double M2_soil;
	
	//LAI
	double mu0_LAI;
	double LAI_MIN;
	double LAI_MAX;
	double T0_LAI;
	double kdecay_LAI;
	double LAI__AT0;
	
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

