/*
 *  spec_model.h
 *  iwq_quality_traditional
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/TRADITIONAL POLLUTANTS
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
	double F_X;
	double C_X;
	double F_wwtp;
	double F_cso;
	double F_direct_storm;
	double F_diffuse;
	double F_diffuse_dissolved;
	double F_diffuse_erosion;
	
	//parameters
	//urban system
	double f_person;
	double C_parasitic;
	double C_storm;
	double K_elim_beta;
	double theta_elim;
	double q_ww;
	
	//rural system
	double C_agro_int;
	double C_agro_ext;
	double C_forest;
	double f_erosion;
	
	//inputs
	double Q_overflow_sewage;
	double Q_overflow_storm;
	double Q_overflow_parasitic;
	double Q_direct_storm;
	double Q_treated_sewage;
	double Q_treated_storm;
	double Q_treated_parasitic;
	double Q_agro_int;
	double Q_agro_ext;
	double Q_forest;
	double T_air;
	double Q_total;
	double rel_tau;
	double F_erosion;
	
public:
	IWQ_MODEL_NAME ();
	virtual void modelFunction(double x); 
	virtual ~IWQ_MODEL_NAME(){ }
	virtual bool verifyParameters();
	virtual bool isStatic(){ return true; }	
};

//------------------------------------------------------------------------------------------

#endif

