/*
 *  spec_model.cpp
 *  iwq_quality_traditional
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/TRADITIONAL POLLUTANTS
 *
 */

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//functions from mathutils.h
double constrain_max(double x, double max);	
double constrain_min(double x, double min);
double constrain_minmax(double x, double min, double max);
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//boundary fluxes
	BFX(F_X);
	BFX(C_X);
	BFX(F_wwtp);
	BFX(F_cso);
	BFX(F_direct_storm);
	BFX(F_diffuse);
	BFX(F_diffuse_erosion);
	BFX(F_diffuse_dissolved);
							
	//parameters
	//urban system
	PAR(f_person);
	PAR(C_parasitic);
	PAR(C_storm);
	PAR(K_elim_beta);
	PAR(theta_elim);
	PAR(q_ww);
	
	//rural system
	PAR(C_agro_int);
	PAR(C_agro_ext);
	PAR(C_forest);
	PAR(f_erosion);
	
	//inputs
	INP(Q_overflow_sewage);
	INP(Q_overflow_storm);
	INP(Q_direct_storm);
	INP(Q_overflow_parasitic);
	INP(Q_treated_sewage);
	INP(Q_treated_storm);
	INP(Q_treated_parasitic);
	INP(Q_agro_int);
	INP(Q_agro_ext);
	INP(Q_forest);
	INP(T_air);
	INP(Q_total);
	INP(rel_tau);
	INP(F_erosion);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	
    double Q_overflow_sewage_eff = constrain_min(Q_overflow_sewage, 0.0);
    double Q_overflow_parasitic_eff = constrain_min(Q_overflow_parasitic, 0.0);
    double Q_overflow_storm_eff = constrain_min(Q_overflow_storm, 0.0);
    double Q_direct_storm_eff = constrain_min(Q_direct_storm, 0.0);
    double Q_treated_parasitic_eff = constrain_min(Q_treated_parasitic, 0.0);
    double Q_treated_sewage_eff = constrain_min(Q_treated_sewage, 0.0);
    double Q_treated_storm_eff = constrain_min(Q_treated_storm, 0.0);
    double Q_agro_int_eff = constrain_min(Q_agro_int, 0.0);
    double Q_agro_ext_eff = constrain_min(Q_agro_ext, 0.0);
    double Q_forest_eff = constrain_min(Q_forest, 0.0);
    double Q_total_eff = 0.0; 
    double q_ww_eff = constrain_min(q_ww, 0.0);
    
    //WASTEWATER SYSTEM
	double rel_tau_comp = rel_tau * pow(theta_elim, T_air - 20.0);
	double K_elim_act = (1.0 - pow(exp(-K_elim_beta), rel_tau_comp)); //was 1.0 - exp(...)
    K_elim_act = constrain_minmax(K_elim_act, 0.0, 1.0);
	
	double C_raw_sewage = constrain_min(q_ww_eff>0.0 ? f_person / q_ww_eff : 0.0, 0.0);
	
	//wwtp
	double F_wwtp_in = Q_treated_sewage_eff * C_raw_sewage + Q_treated_storm_eff * C_storm + Q_treated_parasitic_eff * C_parasitic;
	F_wwtp = (1.0 - K_elim_act) * F_wwtp_in;
	
	//cso
	F_cso = Q_overflow_sewage_eff * C_raw_sewage + Q_overflow_storm_eff * C_storm + Q_overflow_parasitic_eff * C_parasitic;
	
	//direct stormwater
	F_direct_storm = Q_direct_storm_eff * C_storm;
	
	//rural
	F_diffuse_erosion = constrain_min(F_erosion * f_erosion, 0.0);
	F_diffuse_dissolved = constrain_min(Q_agro_int_eff * C_agro_int + Q_agro_ext_eff * C_agro_ext + Q_forest_eff * C_forest, 0.0);
	F_diffuse = F_diffuse_dissolved + F_diffuse_erosion;
	
	F_X = F_wwtp + F_cso + F_diffuse + F_direct_storm;
	
    Q_total_eff = Q_treated_sewage_eff + Q_treated_storm_eff + Q_treated_parasitic_eff + Q_overflow_sewage_eff + Q_overflow_storm_eff +Q_overflow_parasitic_eff + Q_direct_storm_eff + Q_agro_int_eff + Q_agro_ext_eff + Q_forest_eff;
    
	C_X = Q_total_eff>0.0 ? F_X/Q_total_eff : 0.0;
		
	*F(F_X)=F_X;
	*F(C_X)=C_X;
	*F(F_wwtp)=F_wwtp;
	*F(F_cso)=F_cso;
	*F(F_diffuse)=F_diffuse;
	*F(F_diffuse_erosion)=F_diffuse_erosion;
	*F(F_diffuse_dissolved)=F_diffuse_dissolved;
	*F(F_direct_storm)=F_direct_storm;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (K_elim_beta >= 0.0 &&
			theta_elim >= 0.0 && 
			f_person >= 0.0 && 
			q_ww >= 0.0 &&
			C_storm >= 0.0 && 
			C_parasitic >= 0.0 && 
			C_agro_int >= 0.0 && 
			C_agro_ext >=0.0 &&
			C_forest >= 0.0 );
}

//-------------------------------------------------------------------------------------------------------------
