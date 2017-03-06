/*
 *  spec_model.cpp
 *  iwq_heat_soil
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/SOIL
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
	VAR(T_soil);
	VAR(LAI);
	
	//parameters
	//NOTE: endings with "__AT0" will be replaced to "@0" to enable reference to initial values
	PAR(K_soil);
	PAR(M2_soil);
	PAR(mu0_LAI);
	PAR(LAI_MIN);
	PAR(LAI_MAX);
	PAR(T0_LAI);
	PAR(kdecay_LAI);
	PAR(LAI__AT0);
		
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
	double dlai = LAI * (mu0_LAI * constrain_min(T_air - T0_LAI,0) * (1.0 - LAI/LAI_MAX) - kdecay_LAI * constrain_min(T0_LAI - T_air,0));
	if(dlai>0.0){
		*d(LAI) = dlai;
	}
	else{
		double LAI_surplus = LAI - LAI_MIN;
		*d(LAI) =  dlai * (LAI_surplus / (LAI_surplus + 0.1));	//artificial slowdown around LAI_MIN
	}
	double dtsoil = M2_soil * (T_air - T_soil) * ( T_air > T_soil ? exp(-K_soil * LAI) : 1.0);
	if(dtsoil>0.0 || T_soil>0.0){
		*d(T_soil) = dtsoil;
	}
	else{
		*d(T_soil) = 0.0;	//do not allow T_soil to go below 0.0
	}
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	if(LAI_MAX<=0.0 || LAI_MIN>LAI_MAX || mu0_LAI<0.0 || kdecay_LAI<0.0 || M2_soil<0.0 || K_soil<0.0 || LAI__AT0 <0 || LAI__AT0 > LAI_MAX || LAI__AT0 < LAI_MIN){
		return false;
	}
	return true;
}

//-------------------------------------------------------------------------------------------------------------
