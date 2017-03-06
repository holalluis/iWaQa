/*
 *  spec_model.h
 *  iwq_hydrology_terrain
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/TERRAIN
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
	double soil;
	
	//boundary fluxes
	double et;
	double runoff;
	double ssf;
	double rge;
	
	//parameters
	double area;
	double petMult;
	double FC;
	double FS;
	double WP;
	double leachMax;
	double prop_ssf;
	
	//inputs
	double rain_mm;
	double rain_m3s;
	double pet;
	
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
};

//------------------------------------------------------------------------------------------

#endif

