/*
 *  spec_model.cpp
 *  iwq_hydrology_channel
 *
 *  iWaQa model framework 2010-2017
 *
 *  HYDROLOGY/CHANNELS
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//various mathematical utilities (from mathutils.h)
double constrain_min(double x, double min);
double constrain_minmax(double x, double min, double max);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//variables 
	VAR(gw);
	VAR(surf);
	
	//boundary fluxes
	BFX(q);
	BFX(q_new);
	BFX(bf);
	
	//input data
	INP(runoff);
	INP(ssf);
	INP(rge);
	INP(qin);
	
	//parameters
	PAR(kBf);
	PAR(kStream);
	PAR(area);
	PAR(rgeMult);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
    double areaconv = (area / 86.4); 							//mm/d to m3/s if our area unit is km2
	double inverseareaconv = (area!=0.0? 86.4 / area : 0.0);	//m3/s back to mm/d
	double rgeMult_act = constrain_min(rgeMult, 0.0);
	
	double kbf_act = constrain_min(kBf, 0.0);
	double kstream_act = constrain_minmax(kStream, 0.0, 80.0);	
	
	bf = kbf_act * gw;
    q = kstream_act * surf;		//constrained kStream to avoid numerical problems
	
    //derivatives
    *(d(gw)) = rgeMult_act * rge * inverseareaconv - bf;
    *(d(surf)) = qin * inverseareaconv + runoff * inverseareaconv + ssf * inverseareaconv + bf - q;
	
    //fluxes
    *(F(q)) = q * areaconv;
	*(F(q_new)) = q * areaconv - qin;
	*(F(bf)) = bf * areaconv;
}

//-------------------------------------------------------------------------------------------------------------

