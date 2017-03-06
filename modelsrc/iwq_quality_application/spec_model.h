/*
 *  spec_model.h
 *  iwq_quality_application
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/PESTICIDE APPLICATION
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
	double M_to_apply;
	double T_sum;
	
	//boundary fluxes
	double f_applic;
	
	//parameters
	double T_threshold;
	double T_objective;
	double M_total;
	double f_daily;
	double rain_threshold;
	
	//inputs
	double T_air;
	double rain;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return false; }	
};

//------------------------------------------------------------------------------------------

#endif

/*

**************** PESTICIDE APPLICATION ********************

- application switch timer with heat sum (T above threshold accumulates) + total application amount (SUM and PER DAY ALLOWANCE) when condition is met
	- parameters:
		T_threshold;
		T_objective;
		M_total;
		F_daily;
	- variables:
		T_sum;
		M_to_apply;
	- inputs:
		T_air;	
	- outputs:
		F_applic;
- dynamic stock model
	- parameters:
		k_driver; (k_driver * F_driver * M_stock = F_X)
		k_decay; 
	- variables:
		M_stock;
	- inputs:
		F_applic;
		F_driver;
		Q_total;
	- outputs:
		F_X;
		C_X;
		F_decay;
*/

