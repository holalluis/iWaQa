/*
 *  spec_model.cpp
 *  iwq_quality_biocide
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/BIOCIDES
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//functions from mathutils.h
double constrain_max(double x, double max);
double constrain_minmax(double x, double min, double max);
double constrain_min(double x, double min);
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter)
	//variables
	PAR(M_stock);				//stock
	
	//boundary fluxes
	BFX(F_X);				//mobilized flux
	BFX(C_X);				//concentration in streams
	
	//parameters
	PAR(beta);				//washout rate
	PAR(a_urban);
	PAR(C_background);
	PAR(k_wwtp);
		
	//inputs
	INP(F_driver);			//driver flux
	INP(Q_total);			//total stream discharge
	INP(T_air);				//air temperature
	INP(Q_background);		//gw + ssf
	INP(f_wwtp);
	
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	double f_wwtp_eff = constrain_minmax(f_wwtp, 0.0, 1.0);
	F_X = beta * a_urban * F_driver * 0.001 * M_stock * (1.0 - f_wwtp_eff * k_wwtp) + C_background * Q_background * 1E-9; 
	C_X = Q_total>0.0 ? F_X/Q_total : 0.0;
	
	*F(F_X) = F_X;								//kg
	*F(C_X) = C_X * 1E9;						//to µg/m3 = ng/l
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (beta >= 0.0 && k_wwtp >= 0.0 && k_wwtp <= 1.0);
}

//-------------------------------------------------------------------------------------------------------------
