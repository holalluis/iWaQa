/*
 *  spec_model.cpp
 *  iwq_hydrology_canopy
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/CANOPY
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//various mathematical utilities (from mathutils.h)
double SoftMaximum(double x, double y, double k);
double constrain_minmax(double x, double min, double max);
double constrain_min(double x, double min);
double constrain_max(double x, double max);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//variables 
	VAR(intercept_storage);
	
	//boundary fluxes
	BFX(et_mm);
	BFX(throughfall_mm);
	BFX(et_m3s);
	BFX(throughfall_m3s);
	
	//input data
	INP(rain_mm);
	INP(pet_mm);
	
	//parameters
	PAR(harvest_eff);
	PAR(storage_size);
	PAR(LAI_min);
	PAR(petMult);
	PAR(area);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	double areaconv = (area / 86.4); 							//mm/d to m3/s if our area unit is km2
	
	//parameter filtering
	double harvest_eff_act = constrain_minmax(harvest_eff, 0.0, 1.0);
	double storage_size_act = constrain_min(storage_size, 0.0);
	double LAI_min_act = constrain_minmax(LAI_min, 0.0, 1.0);
	double petMult_act = constrain_min(petMult, 0.0);
	
	//simplified LAI substitute
	double DOY = x - (365.0 * int( x / 365.0));
	double sine = sin(DOY * M_PI / 365.0);
	double LAI_act = LAI_min_act + (1.0 - LAI_min_act) * sine * sine;
	
	//flux calculations  
	double harvested = rain_mm * harvest_eff_act * LAI_act;
	double threshold_storage = storage_size_act * LAI_act;
	double leaked = 86.4 * constrain_min( intercept_storage - threshold_storage, 0.0); //SoftMaximum( intercept_storage - threshold_storage, 0.0, 50.0);	//hardcoded 1000s residence time in canopy
	double evaporated = pet_mm * petMult_act * LAI_act * (intercept_storage / (intercept_storage + 0.1));
	
	//underflow guard
	if(intercept_storage<0.0){
		leaked=0.0;
		evaporated=0.0;
	}
	
	//derivatives
	*(d(intercept_storage)) = harvested - leaked - evaporated;
	
	double tf = (leaked + (rain_mm - harvested));	//throughfall
	
    //fluxes (converted to m3/s here, above everything is in mm/d)
    *(F(et_mm)) = evaporated;
    *(F(throughfall_mm)) = tf;	//including immediate throughfall
	*(F(et_m3s)) = evaporated * areaconv;
    *(F(throughfall_m3s)) = tf * areaconv;	//including immediate throughfall
}

//-------------------------------------------------------------------------------------------------------------

