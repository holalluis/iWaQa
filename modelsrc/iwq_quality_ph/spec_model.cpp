/*
 *  spec_model.cpp
 *  iwq_quality_ph
 *
 *  iWaQa model framework 2010-2017
 *
 *  QUALITY/pH
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
double SoftThreshold(double x, double threshold, double k);
double SoftMaximum(double x, double y, double k);

//############################################################################################################

IWQ_MODEL_NAME::IWQ_MODEL_NAME() : iWQModel( QUOTEME( IWQ_MODEL_NAME ) )
{	
	// NOTE: specify types with VAR (variable), BFX (boundary flux), INP (input) and PAR (parameter) 
	//boundary fluxes
	BFX(pH);
	BFX(C_TIC);
	BFX(F_alk);
	BFX(C_alk);
			
	//parameters
	PAR(pH_wwtp);
	PAR(pH_raw_sewage);
	PAR(pH_natural);
	PAR(pH_rain);
	PAR(C_alk_wwtp);
	PAR(C_alk_raw_sewage);
	PAR(C_alk_natural);
			
	//inputs
	INP(Q_overflow_sewage);
	INP(Q_overflow_storm);
	INP(Q_overflow_parasitic);
	INP(Q_wwtp);
	INP(T_air);
	INP(Q_total);

}

//-------------------------------------------------------------------------------------------------------------

double pKa(double temp)
{
	double pka = 6.57 - 0.0118 * temp + 0.00012 * (temp * temp);
	return pka;
}

double ionfrac(double pH, double temp)
{
	return 1.0 / (1.0 + pow(10.0, pKa(temp) - pH));
}

//-------------------------------------------------------------------------------------------------------------

void IWQ_MODEL_NAME::modelFunction(double x)
{
	// NOTE: calculate changes, then assign with *(d()) (derivative, for VARs) and *(F()) (flux, for BFXs)
	pH=7.0;
	F_alk=0.0;
	C_alk=0.0;
	C_TIC=0.0;
	
	if(Q_total>0.0){
			
		//MIX CSO EFFLUENT COMPONENTS
		double Q_CSO = (Q_overflow_sewage + Q_overflow_storm + Q_overflow_parasitic);
		double C_alk_CSO = 0.0;
		double pH_CSO = 7.0;
		double TIC_CSO = 0.0;
		
		if(Q_CSO>0.0){
			double TIC_overflow_parasitic = C_alk_natural / ionfrac(pH_natural, T_air);
			double TIC_overflow_sewage = C_alk_raw_sewage / ionfrac(pH_raw_sewage, T_air);
		
			C_alk_CSO = (Q_overflow_sewage * C_alk_raw_sewage + Q_overflow_parasitic * C_alk_natural) / Q_CSO;
			TIC_CSO = (Q_overflow_sewage * TIC_overflow_sewage + Q_overflow_parasitic * TIC_overflow_parasitic) / Q_CSO;
			pH_CSO = pKa(T_air) - log10((TIC_CSO/C_alk_CSO)-1.0);
		}
		
		//MIX WWTP EFFLUENT WITH THE NATURAL FLOW (MIX1)
		double Q_nat = Q_total - Q_CSO - Q_wwtp;
		
		double TIC_nat = C_alk_natural / ionfrac(pH_natural, T_air);
		double TIC_wwtp = C_alk_wwtp / ionfrac(pH_wwtp, T_air);
		
		double TIC_mix1 = 0.0;
		double pH_mix1 = 7.0;
		double C_alk_mix1 = 0.0;
		double Q_mix1 = Q_nat + Q_wwtp;
		if(Q_mix1 > 0.0){
			TIC_mix1 = (Q_nat * TIC_nat + Q_wwtp * TIC_wwtp) / Q_mix1;
			C_alk_mix1 = (Q_nat * C_alk_natural + Q_wwtp * C_alk_wwtp) / Q_mix1;
			pH_mix1 = pKa(T_air) - log10((TIC_mix1/C_alk_mix1)-1.0);
		}
		
		//MIX MIX1 WITH CSO
		C_TIC = (Q_mix1 * TIC_mix1 + Q_CSO * TIC_CSO) / Q_total;
		C_alk = (Q_mix1 * C_alk_mix1 + Q_CSO * C_alk_CSO) / Q_total;
		pH = pKa(T_air) - log10((C_TIC/C_alk)-1.0);
	}
	
	*F(pH) = pH;
	*F(C_TIC) = C_TIC;
	*F(C_alk) = C_alk;
	*F(F_alk) = F_alk;
}

//-------------------------------------------------------------------------------------------------------------

bool IWQ_MODEL_NAME::verifyParameters()
{
	//return here false if the actual parameter configuration contains illegal values
	return (pH_wwtp >= 0.0 && pH_wwtp <= 14.0 &&
			pH_raw_sewage >= 0.0 && pH_raw_sewage <= 14.0 &&
			pH_natural >= 0.0 && pH_natural <= 14.0 &&
			pH_rain >= 0.0 && pH_rain <= 14.0 &&
			C_alk_wwtp >= 0.0 && 
			C_alk_raw_sewage >= 0.0 && 
			C_alk_natural >= 0.0);
}

//-------------------------------------------------------------------------------------------------------------
