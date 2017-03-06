/*
 *  spec_model.h
 *  iwq_quality_biocide
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/BIOCIDES
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
	double M_stock; 		//stock
	
	//boundary fluxes
	double F_X;				//mobilized flux
	double C_X;				//concentration in streams
	
	//parameters
	double beta;			//washout rate
	double a_urban;
	double C_background;
	double k_wwtp;			//constant elimination rate at WWTP
		
	//inputs
	double F_driver;		//driver flux
	double f_wwtp;			//proportion of F_driver passing through the WWTP
	double Q_total;			//total stream discharge
	double T_air;			//air temperature
	double Q_background;	//subsurface water flux (gw + ssf)
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return true; }	//<-- STATIC MODEL!
};

//------------------------------------------------------------------------------------------

#endif

