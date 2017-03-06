/*
 *  spec_model.h
 *  iwq_heat_mixer
 *
 *  iWaQa model framework 2010-2017
 *
 *  HEAT/MIXER
 *
 */

#ifndef spec_model_h
#define spec_model_h

#include "model.h"

//--------------------------------------------------------------------------------------------

class IWQ_MODEL_NAME : public iWQGenericChannelTransport
{
private:
	// NOTE: declare here every variable, input, flux and parameter as double

	//variables
	double T_water;
	
public:
	IWQ_MODEL_NAME ();
	virtual ~IWQ_MODEL_NAME(){ }
	virtual double R(double x);
	virtual bool verifyParameters();
};

//------------------------------------------------------------------------------------------

#endif

