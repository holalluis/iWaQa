/*
 *  spec_model.cpp
 *  iwq_heat_network
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/STREAM NETWORK
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	BFX(T_water);
	BFX(T_water_flux);
	
	//parameters
	PAR(tauperd);
	PAR(T_offset);	//for compensation of altitude differences. Assumption: T_source and T_eq are affected by the same offset due to a difference in T_air 
	PAR(shaded);
	PAR(LAI_MAX);
		
	//inputs
	INP(T_air);
	INP(T_source);
	INP(T_eq);
	INP(T_eq_shade);
	INP(K_model);
	INP(K_model_shade);
	INP(Q);
	INP(LAI);
}

//-------------------------------------------------------------------------------------------------------------

//handy functions from mathutils.h
double constrain_minmax(double x, double min, double max);	//equals to min(max,max(x,min))
double constrain_min(double x, double min);					//equals to max(x,min)
double constrain_max(double x, double max);					//equals to min(x,max)

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	double shaded_eff = shaded * constrain_max(LAI / LAI_MAX, 1.0);
	double T_eq_eff = shaded_eff * T_eq_shade + (1.0 - shaded_eff) * T_eq; 
	double K_model_eff = shaded_eff * K_model_shade + (1.0 - shaded_eff) * K_model; 
	T_water = T_eq_eff + (T_source - T_eq_eff) * exp( -K_model_eff * tauperd ) + T_offset;
	*F(T_water) = T_water;
	*F(T_water_flux)= Q * T_water; //heat flux
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	if(tauperd<0.0 || K_model<0.0 || K_model_shade<0.0 || shaded<0.0 || shaded>1.0 || LAI_MAX<=0.0){
		return false;
	}
	return true;
}

//-------------------------------------------------------------------------------------------------------------
