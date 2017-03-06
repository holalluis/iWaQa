/*
 *  spec_model.h
 *  iwq_quality_pesticide
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/PESTICIDE TRANSPORT
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
	double M_stock;
	
	//boundary fluxes
	double F_X;
	double C_X;
	double F_decay;
	
	//parameters
	double beta;
	double appl_loss;
	double k_decay;
	double theta_decay;
	double area;
	double area_applic;
	double C_background;
	
	//inputs
	double F_driver;
	double f_applic;	//application flux
	double Q_total;
	double T_air;
	double Q_background;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return false; }	
};

//------------------------------------------------------------------------------------------

#endif

