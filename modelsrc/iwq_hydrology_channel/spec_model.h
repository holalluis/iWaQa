/*
 *  spec_model.h
 *  iwq_hydrology_channel
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/CHANNELS
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
	double gw;
	double surf;
	
	//boundary fluxes
	double q;
	double q_new;
	double bf;
	
	//parameters
	double area;
	double kBf;
	double kStream;
	double rgeMult;
	
	//inputs
	double runoff;
	double ssf;
	double rge;
	double qin;

public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
};

//------------------------------------------------------------------------------------------

#endif

