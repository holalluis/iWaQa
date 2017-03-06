/*
 *  spec_model.cpp
 *  iwq_heat_mixer
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/MIXER
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//***********************************************************************************************************

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQGenericChannelTransport( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with CVAR (concentration variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//conc. variable
	CVAR(T_water);
}

//-------------------------------------------------------------------------------------------------------------

double IWQ_MODEL_NAME::R(double x)		
{
	//return the reaction rate for the concentration variable
	return 0.0;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return true;
}

//-------------------------------------------------------------------------------------------------------------


