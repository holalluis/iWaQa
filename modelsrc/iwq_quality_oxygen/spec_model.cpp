/*
 *  spec_model.cpp
 *  iwq_quality_oxygen
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/DISSOLVED OXYGEN
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
	//boundary fluxes
	BFX(F_DO);
	BFX(C_DO);
	BFX(C_DO_sat);
			
	//parameters
	PAR(C_wwtp);
	PAR(C_raw_sewage);
	PAR(C_parasitic);
	PAR(elevation);
	PAR(T_water_min);
			
	//inputs
	INP(Q_overflow_sewage);
	INP(Q_overflow_storm);
	INP(Q_overflow_parasitic);
	INP(Q_wwtp);
	INP(T_air);
	INP(Q_total);

}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
	double T_water = constrain_min(T_air, T_water_min);
	
	C_DO_sat = 14.6 * exp(-(0.027767+(-0.00027+0.000002*T_water)*T_water)*T_water)*pow(1.0-0.00000697*elevation*3.28,5.167);
	
	double Q_cso = Q_overflow_sewage + Q_overflow_storm + Q_overflow_parasitic;
	double Q_natural = Q_total - Q_wwtp - Q_cso;
	
	double F_cso = Q_overflow_sewage * C_raw_sewage + Q_overflow_storm * C_DO_sat + Q_overflow_parasitic * C_parasitic;
	double F_wwtp = Q_wwtp * C_wwtp;
	double F_natural = Q_natural * C_DO_sat;
	
	F_DO = F_wwtp + F_cso + F_natural;
	C_DO = Q_total>0.0 ? F_DO/Q_total : 0.0;
		
	*F(F_DO)=F_DO;
	*F(C_DO)=C_DO;
	*F(C_DO_sat)=C_DO_sat;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (C_raw_sewage >= 0.0 && 
			C_parasitic >= 0.0 && 
			C_wwtp >= 0.0);
}

//-------------------------------------------------------------------------------------------------------------
