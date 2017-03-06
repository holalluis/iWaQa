/*
 *  spec_model.cpp
 *  iwq_heat_wastewater
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/WASTEWATER
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
	VAR(T_air_smooth);
	BFX(T_ww);
	
	//parameters
	PAR(W_smoothing);
	PAR(T_source_ww);
	PAR(K_ww);
		
	//inputs
	INP(T_air);
}

//-------------------------------------------------------------------------------------------------------------

//handy functions from mathutils.h
double constrain_minmax(double x, double min, double max);	//equals to min(max,max(x,min))
double constrain_min(double x, double min);					//equals to max(x,min)
double constrain_max(double x, double max);					//equals to min(x,max)

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	*d(T_air_smooth) = 1.0 / W_smoothing * (T_air - T_air_smooth);
	*F(T_ww) = T_air_smooth + (T_source_ww - T_air_smooth) * K_ww;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	if(W_smoothing<=0.0 || K_ww<0.0){
		return false;
	}
	return true;
}

//-------------------------------------------------------------------------------------------------------------
