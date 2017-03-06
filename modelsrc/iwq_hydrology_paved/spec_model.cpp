/*
 *  spec_model.cpp
 *  iwq_hydrology_paved
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/PAVED SURFACES
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

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//variables 
	VAR(storage);
	
	//boundary fluxes
	BFX(runoff);
	BFX(et);
	BFX(infiltration);
	
	//input data
	INP(rain);
	INP(pet);
	
	//parameters
	PAR(area);
	PAR(s);
    PAR(s_mult);
	PAR(k_s);
	PAR(petMult);
	PAR(k_infiltr);
    PAR(k_impermeable);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
    double areaconv = (area / 86.4); 							//mm/d to m3/s if our area unit is km2
	double inverseareaconv = (area!=0.0? 86.4 / area : 0.0);	//m3/s back to mm/d
	
	double s_eff=constrain_min(s * s_mult, 0.0);
	double k_s_eff=constrain_minmax(k_s, 0.0, 20.0);
	
	runoff=SoftMaximum(k_s_eff * (storage-s_eff), 0.0, 5.0);
	if(runoff>storage){
		runoff=storage;
	}
	
	et = (petMult * pet < storage) ? petMult * pet * (storage / (storage+0.1)) : storage;	//mm	
	
	infiltration = (1.0 - k_impermeable) * k_infiltr * storage / (storage + 0.1);
	
    //derivatives
    *(d(storage)) = rain - runoff - et - infiltration;
    
    //fluxes
    *F(et) = et * areaconv;		//fluxes are converted to m3/s here, above everything is in mm/d
    *F(runoff) = runoff * areaconv;
	*F(infiltration) = infiltration * areaconv;
}

//-------------------------------------------------------------------------------------------------------------

