/*
 *  spec_model.cpp
 *  iwq_quality_pesticide
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/PESTICIDE TRANSPORT
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//functions from mathutils.h
double constrain_max(double x, double max);	
double constrain_min(double x, double min);
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter)
	//variables
	VAR(M_stock);
	
	//boundary fluxes
	BFX(F_X);
	BFX(C_X);
	BFX(F_decay);
	
	//parameters
	PAR(beta);	//
	PAR(appl_loss);
	PAR(k_decay);
	PAR(theta_decay);
	PAR(area);
	PAR(area_applic);	//proportion of area where the compound is applied
	PAR(C_background);
		
	//inputs
	INP(F_driver);	//driver hydrologic flux (m3)
	INP(f_applic);	//application flux per area unit (km2)
	INP(Q_total);
	INP(T_air);
	INP(Q_background);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	F_decay = k_decay * M_stock * pow(theta_decay, T_air - 20.0);
	double F_X_stock = beta * F_driver/area * 0.001 * M_stock + f_applic*area_applic*appl_loss;	// [prop/mm] * [m3]/[km2] * 0.001[mm/(m3/km2)] * kg
	double F_X_bg = Q_background * C_background * 1e-9;
	F_X = F_X_stock + F_X_bg;
	C_X = Q_total>0.0 ? F_X/Q_total : 0.0;
	
	*d(M_stock) = f_applic * area_applic - F_decay - F_X_stock;	//f_applic is in kg/km2 while we want kg
	*F(F_decay) = F_decay;	//kg
	*F(F_X) = F_X;			//kg
	*F(C_X) = C_X * 1E9;	//to µg/m3 = ng/l
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (k_decay >= 0.0 && beta >= 0.0 && theta_decay >= 0.5 && theta_decay < 1.5 && appl_loss>=0.0);
}

//-------------------------------------------------------------------------------------------------------------
