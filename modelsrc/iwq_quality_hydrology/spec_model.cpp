/*
 *  spec_model.cpp
 *  iwq_quality_hydrology
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/POLLUTANT HYDROLOGY
 *
 */
 
//HYDROLOGICAL FLUXES FOR WATER QUALITY

#include "model.h"

#include "spec_model.h"

#include <math.h>

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

//functions from mathutils.h
double constrain_max(double x, double max);	
double constrain_min(double x, double min);
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	
	//boundary fluxes
	BFX(Q_wwtp);
	BFX(Q_cso);
	BFX(Q_storm);
	BFX(Q_storm_direct);
	BFX(rel_tau);
	BFX(Q_overflow_sewage);
	BFX(Q_overflow_storm);
	BFX(Q_overflow_parasitic);
	BFX(Q_treated_sewage);
	BFX(Q_treated_storm);
	BFX(Q_treated_parasitic);
	BFX(q_diffuse);
	BFX(Q_agro_int);
	BFX(Q_agro_ext);
	BFX(Q_forest);
	BFX(F_erosion);
	
	//parameters
	PAR(k_storm_runoff);	//proportion of Q_runoff coming in the sewers
	PAR(k_storm_direct_runoff); //proportion of Q_runoff coming in the separated sewers
	PAR(k_storm_ssf);		//proportion of Q_ssf coming in the sewers
	PAR(Q_cso_threshold);	//maximal incoming Q_storm to WWTP
	PAR(k_cso);				//shape parameter for the overflow curve
	PAR(q_ww);
	PAR(n_person);
	PAR(k_mixed_flow);
	PAR(A_total);
	PAR(A_agro_int);
	PAR(A_agro_ext);
	PAR(A_forest);
	PAR(a_erosion);
	PAR(b_erosion);
	
	//inputs
	INP(Q_runoff);
	INP(Q_ssf);
	INP(Q_total);
	INP(rain);
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	//derivatives
	
	//WASTEWATER SYSTEM
	
	//calculate diffuse and point source discharges
	Q_wwtp = 0.0;			//WWTP hydraulic load
	Q_cso = 0.0;			//CSO discharge
	
	//domestic wastewater
	double Q_domestic = n_person * q_ww;
	
	//parasitic water
	double Q_parasitic = k_storm_ssf * Q_ssf;
	
	//steady sewer flow
	double Q_steady = Q_domestic + Q_parasitic;
	
	//default non-storm case
	Q_treated_sewage = Q_domestic;
	Q_treated_parasitic = Q_parasitic;
	Q_treated_storm = 0.0;
	Q_overflow_sewage = 0.0;
	Q_overflow_storm = 0.0;
	Q_overflow_parasitic = 0.0;
	
	//storm runoff & overflow (mixed system)
	Q_storm = k_storm_runoff * Q_runoff;
	if(Q_storm>0.0){
		double Q_total_combined_sewer = k_mixed_flow * Q_steady + Q_storm;
		double Q_intake = -SoftMaximum(-Q_total_combined_sewer, -Q_cso_threshold, k_cso); 
		Q_cso = Q_total_combined_sewer - Q_intake;
		double p_overflow = Q_cso / Q_total_combined_sewer;
		Q_overflow_sewage = p_overflow * k_mixed_flow * Q_domestic;
		Q_overflow_storm = p_overflow * Q_storm;
		Q_overflow_parasitic = p_overflow * k_mixed_flow * Q_parasitic;
		Q_treated_sewage = Q_domestic - Q_overflow_sewage;
		Q_treated_storm = Q_storm - Q_overflow_storm;
		Q_treated_parasitic = Q_parasitic - Q_overflow_parasitic;
	}
	
	Q_wwtp = Q_steady + Q_storm - Q_cso;
	
	*F(Q_wwtp) = Q_wwtp;							// [m3 d-1]
	*F(Q_cso) = Q_cso;								// [m3 d-1]
	*F(Q_storm) = Q_storm;							// [m3 d-1]
	*F(Q_overflow_sewage)=Q_overflow_sewage;		// [m3 d-1]
	*F(Q_overflow_storm)=Q_overflow_storm;			// [m3 d-1]
	*F(Q_overflow_parasitic)=Q_overflow_parasitic;	// [m3 d-1]
	*F(Q_treated_sewage)=Q_treated_sewage;			// [m3 d-1]
	*F(Q_treated_storm)=Q_treated_storm;			// [m3 d-1]
	*F(Q_treated_parasitic)=Q_treated_parasitic;	// [m3 d-1]
	
	//relative residence time
	rel_tau = Q_wwtp > 0.0 ? Q_domestic / Q_wwtp : 1.0; 
	
	*F(rel_tau) = rel_tau;	// [d d-1]
	
	//SEPARATE STORM SEWER SYSTEM
	Q_storm_direct = k_storm_direct_runoff * Q_runoff;
	
	*F(Q_storm_direct) = Q_storm_direct;			// [m3 d-1]
	
	//RURAL (DIFFUSE) SYSTEM
	
	//fluxes for soluble pollutants
	double Q_rural = constrain_min(Q_total - Q_storm - Q_steady - Q_storm_direct, 0.0);
	double A_rural = A_agro_int + A_agro_ext + A_forest;
	q_diffuse = A_rural>0.0 ? Q_rural / A_rural : 0.0;
	Q_agro_int = A_agro_int * q_diffuse;
	Q_agro_ext = A_agro_ext * q_diffuse;
	Q_forest = A_forest * q_diffuse;
	
	*F(q_diffuse)=q_diffuse;
	*F(Q_agro_int)=Q_agro_int;
	*F(Q_agro_ext)=Q_agro_ext;
	*F(Q_forest)=Q_forest;
	
	//erosion flux (virtual material flux from intensive agriculture)
	F_erosion = A_agro_int * a_erosion * pow(rain, b_erosion);
	
	*F(F_erosion)=F_erosion;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (k_storm_runoff >= 0.0 && k_storm_runoff <= 1.0 &&
			k_storm_ssf >= 0.0 && k_storm_ssf <= 1.0 &&
			k_mixed_flow >= 0.0 && k_mixed_flow <= 1.0 &&
			Q_cso_threshold >= 0.0 && k_cso >= 0.0 && q_ww >= 0.0 && n_person>=0.0 &&
			A_total >= 0.0 && A_agro_int >= 0.0 && A_agro_ext >= 0.0 && A_forest >= 0.0 && 
			A_total >= (A_agro_int + A_agro_ext + A_forest));
}

//-------------------------------------------------------------------------------------------------------------
