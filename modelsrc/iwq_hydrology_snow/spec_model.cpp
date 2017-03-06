/*
 *  spec_model.cpp
 *  iwq_hydrology_snow
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/SNOW
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//various mathematical utilities (from mathutils.h)
double SoftMaximum(double x, double y, double k);
double SoftThreshold(double x, double threshold, double k);
double constrain_minmax(double x, double min, double max);
double constrain_min(double x, double min);
double constrain_max(double x, double max);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//variables 
	VAR(snow);
	
	//boundary fluxes
	BFX(rain);
	
	//input data
	INP(prec);
	INP(temp);
	
	//parameters
	PAR(rMult);
	PAR(tcrit);
	PAR(tmelt);
	PAR(ksnow);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
	//parameter filtering
	double rMult_act = constrain_min(rMult, 0.0);
	double ksnow_act = constrain_min(ksnow, 0.0);
	
	double effluent = 0.0;
	double new_snow = 0.0;
	double melt = 0.0;
		
	//Soft rain threshold	
	double rain_proportion = SoftThreshold(temp, tcrit, 1.0);
	
	effluent = rMult_act * prec * rain_proportion;
	new_snow = rMult_act * prec * (1.0 - rain_proportion);
	
	//snowmelt with soft maximum
	if(temp>tmelt){
		melt=SoftMaximum(ksnow_act*(temp-tmelt), 0.0, 5.0);
		melt=(melt>snow?snow:melt);	//limit melt (hard min function)
	}

	//change in snow depth
	*(d(snow)) = new_snow - melt;
    
    //fluxes (here in mm)
    *(F(rain)) = effluent + melt;
}

//-------------------------------------------------------------------------------------------------------------

