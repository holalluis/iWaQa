/*
 *  spec_model.cpp
 *  iwq_hydrology_terrain
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/TERRAIN
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//functions from mathutils.h
double constrain_max(double x, double max);	
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//variables 
	VAR(soil);
	
	//boundary fluxes
	BFX(et);
	BFX(runoff);
	BFX(ssf);
	BFX(rge);
	
	//input data
	INP(rain_mm);
	INP(rain_m3s);
	INP(pet);
	
	//parameters
	PAR(area);
	PAR(petMult);
	PAR(FC);
	PAR(FS);
	PAR(WP);
	PAR(leachMax);
	PAR(prop_ssf);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
    double areaconv = (area / 86.4); 							//mm/d to m3/s if our area unit is km2
	double inverseareaconv = (area!=0.0? 86.4 / area : 0.0);	//m3/s back to mm/d
	
	//unify inputs
	double rain = rain_mm + inverseareaconv * rain_m3s;
	
	//notation
	// WP: wilting point
	// FC: field capacity
	// FS: full saturation
	
	//saturated area
	double h_s50 = (FS + FC) / 2.0;
	double sigma= (FS - FC) / 4.0;
	
	double f_sat = (1.0 / (1.0 + exp(2.0/sigma * (h_s50 - soil))) - 1.0 / (1.0 + exp( 2 * h_s50 / sigma))) ;	//saturated area function
	    
	runoff = rain * f_sat;
	
	ssf = leachMax * prop_ssf * f_sat;
    rge = leachMax * (1.0 - prop_ssf) * f_sat;
    
	//evapotranspiration
	double et_50 = 0.25 * (3.0 * WP + FC);		// middle point between WP and (WP+FC)/2, which corresponds to the easily available water content limit (pF=3.3)
	double k_shape = 10.0 / et_50;				// fixed shape factor so that the change happens in the given interval
	double et_per_pet = SoftThreshold(soil, et_50, k_shape) - SoftThreshold(0, et_50, k_shape); 	//assumes that the _limit_ for ET is at (WP+FC)/2
	
	et = pet *  et_per_pet;
	    
    //derivatives
    *(d(soil)) = rain - runoff - ssf - rge - et;
	
    //fluxes (converted to m3/s here, above everything is in mm/d)
    *(F(et)) = et * areaconv;
    *(F(runoff)) = runoff * areaconv;
	*(F(ssf)) = ssf * areaconv;
	*(F(rge)) = rge * areaconv;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (FS > FC && prop_ssf>=0.0 && prop_ssf <=1.0 && area > 0.0 && leachMax >= 0.0);
}
