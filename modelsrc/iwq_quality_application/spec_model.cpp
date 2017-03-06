/*
 *  spec_model.cpp
 *  iwq_quality_application
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/PESTICIDE APPLICATION
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
	VAR(M_to_apply);
	VAR(T_sum);
	
	//boundary fluxes
	BFX(f_applic);
	
	//parameters
	PAR(T_threshold);
	PAR(T_objective);
	PAR(M_total);	// kg / km2
	PAR(f_daily);	// kg / km2
	PAR(rain_threshold);
	
	//inputs
	INP(T_air);
	INP(rain);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
	double t_reset_stock = -5.0;
	double t_reset_tsum = 0.0;
    
    double doy = fmod(x, 365.25);
	
	//heat sum model
	double increment = constrain_min(T_air - T_threshold, 0.0);
	double freeze = (T_air <= t_reset_tsum)? 0.99 * T_sum : 0.0;	//resets the heat sum when it freezes
	*d(T_sum) = increment - freeze;
	
	//application model
	f_applic = 0.0;
    double replenish = 0.0;
	if(doy>=5.0){ 
        if(M_total > 0.0 && T_sum > T_objective && rain<=rain_threshold){
            f_applic = f_daily * M_to_apply / (M_to_apply + 0.01 * M_total);
        }
        *d(M_to_apply) = - f_applic;
    }
    else{
        *d(M_to_apply) = 2.0 * (M_total - M_to_apply);      //resets the application stock (new year)
    }
		
	*F(f_applic) = f_applic;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (T_objective >= 0.0 && M_total >= 0.0 && f_daily >= 0.0 && rain_threshold>=0.0);
}

//-------------------------------------------------------------------------------------------------------------
