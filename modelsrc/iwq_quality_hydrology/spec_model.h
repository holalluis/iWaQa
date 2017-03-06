/*
 *  spec_model.h
 *  iwq_quality_hydrology
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/POLLUTANT HYDROLOGY
 *
 */

#ifndef spec_model_h
#define spec_model_h

#include "model.h"

class iWQModel;

//------------------------------------------------------------------------------------------

class IWQ_MODEL_NAME : public iWQModel
{
private:
	// NOTE: declare here every variable, input, flux and parameter as double

	//boundary fluxes
	double Q_wwtp;
	double Q_cso;
	double Q_storm;
	double Q_storm_direct;
	double Q_overflow_sewage;
	double Q_overflow_storm;
	double Q_overflow_parasitic;
	double Q_treated_sewage;
	double Q_treated_storm;
	double Q_treated_parasitic;
	double rel_tau;
	double q_diffuse;
	double Q_agro_int;
	double Q_agro_ext;
	double Q_forest;
	double F_erosion;
					
	//parameters
	double k_storm_runoff;			//proportion of Q_runoff coming in the mixed sewer system
	double k_storm_direct_runoff; 	//proportion of Q_runoff coming in the separate sewers
	double k_storm_ssf;				//proportion of Q_ssf coming in the mixed sewers (separate sewers are shallow & short)
	double Q_cso_threshold;			//maximal incoming Q_storm
	double k_cso;					//shape parameter for the overflow curve
	double q_ww;
	double n_person;
	double k_mixed_flow;			//proportion of flow mixing (combined sewer proportion from steady Q)
	double A_total;
	double A_agro_int;
	double A_agro_ext;
	double A_forest;
	double a_erosion;
	double b_erosion;

	//inputs
	double Q_runoff;
	double Q_ssf;
	double Q_total;
	double rain;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return true; }
};

//------------------------------------------------------------------------------------------

#endif

