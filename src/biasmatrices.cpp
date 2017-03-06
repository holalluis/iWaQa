/*
 *  biasmatrices.cpp
 *  Helper functions for complicated likelihood calculations
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */

#include <math.h>
#include <iostream> 		//Eigen needs it
#include "biasmatrices.h"
#include "Eigen/Dense"
#include "Eigen/Cholesky"
#include "mathutils.h"

//---------------------------------------------------------------------------------------------------

#pragma mark Bias-specific routines

//RETURNS THE STDEV ADDITION DUE TO INPUT
double g(double input, double kappa, double pi)
{
	if(pi==0.0){
		return input * kappa; 
	}
	else{
		return  pow(input, pi+1.0) * kappa;
	}
}

//---------------------------------------------------------------------------------------------------

//RETURNS THE JUMP VARIANCE OF B
double jumpVarianceOfB(double sigma_b2, double beta, double kappa, double pi, double input)
{
	double inpdep = g(kappa, input, pi);
	return sigma_b2 * (1.0 - exp(-2.0 * beta)) + inpdep * inpdep;
}

//---------------------------------------------------------------------------------------------------

//RETURNS THE VARIANCE OF E
double varianceOfE(double input, double sigma_e2, double kappa_e)
{
	double inpdep = input * kappa_e;
	return sigma_e2 + inpdep * inpdep;
}

//---------------------------------------------------------------------------------------------------

//PREPARES SIGMA_B
Eigen::MatrixXd makeSigmaBMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi)
{
	int md=inputs.size();
	
	Eigen::MatrixXd SIGMA_B = Eigen::MatrixXd::Zero(md, md);
	
	double unconditional_variance=sigma_b2;
	
	//diagonal-controlled filling
	for(int d=0; d<md; d++){
		double jump_var = jumpVarianceOfB(sigma_b2, beta, kappa, pi, inputs[d]);
		unconditional_variance = exp(-2.0 * beta) * unconditional_variance + jump_var;
		SIGMA_B(d,d) = unconditional_variance;
		for(int r=d+1; r<md; r++){
			SIGMA_B(r,d)=SIGMA_B(r-1,d) * exp(-beta);
		}
		for(int c=d+1; c<md; c++){
			SIGMA_B(d,c)=SIGMA_B(d, c-1) * exp(-beta);
		}
	}
	
	return SIGMA_B;
}

//---------------------------------------------------------------------------------------------------

//PREPARES SIGMA_E

Eigen::MatrixXd makeSigmaEMatrix(std::vector<double> inputs, double sigma_e2, double kappa_e)
{
	int md=inputs.size();
	
	Eigen::MatrixXd SIGMA_E = Eigen::MatrixXd::Zero(md, md);
	
	//diagonal filling
	for(int d=0; d<md; d++){
		SIGMA_E(d,d) = varianceOfE(inputs[d], sigma_e2, kappa_e);
	}
	
	return SIGMA_E;
}

//---------------------------------------------------------------------------------------------------

//SUPPLIES SIGMA_E^-1

Eigen::MatrixXd makeSigmaEInverse(std::vector<double> inputs, double sigma_e2, double kappa_e)
{
	int md=inputs.size();
	
	Eigen::MatrixXd SIGMA_EINV = Eigen::MatrixXd::Zero(md, md);
	
	//diagonal filling
	for(int d=0; d<md; d++){
		SIGMA_EINV(d,d) = 1.0 / varianceOfE(inputs[d], sigma_e2, kappa_e);
	}
	
	return SIGMA_EINV;
}

//---------------------------------------------------------------------------------------------------

//INVERTS ANY TRIDIAGONAL MATRIX
Eigen::MatrixXd generalInvertTridiagonal(Eigen::MatrixXd & T)
{
	//Usmani's original method without assuming symmetry
	int n=T.cols();
	Eigen::VectorXd theta (n+1);	//1-based indexing for both
	Eigen::VectorXd phi (n+2);
	
	theta[0]=1.0;
	theta[1]=T(0,0);
	for(int i=2; i<=n; i++){
		theta[i] = T(i-1, i-1) * theta[i-1] - T(i-2, i-1) * T(i-1, i-2) * theta[i-2];
	}
	
	phi[n+1]=1.0;
	phi[n]=T(n-1, n-1);
	for(int i=n-1; i>=1; i--){
		phi[i] = T(i-1, i-1) *  phi[i+1] - T(i-1, i) * T(i, i-1) * phi[i+2];
	}
	
	Eigen::MatrixXd Tinv (n, n);
	//construct lower triangle
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			double minuspart = ((i+j)%2)?-1.0:1.0;
			double prodpart = 1.0;
			
			if(i<=j){
				for(int k=i; k<=j-1; k++){
					prodpart *= T(k-1, k);
				}
				Tinv(i-1, j-1) = minuspart * prodpart * theta[i-1] * phi[j+1] / theta[n];
			}
			else{
				for(int k=j; k<=i-1; k++){
					prodpart *= T(k, k-1);
				}
				Tinv(i-1, j-1) = minuspart * prodpart * theta[j-1] * phi[i+1] / theta[n];
			}
		}
	}
	
	return Tinv;
}

//---------------------------------------------------------------------------------------------------

//SUPPLIES SIGMA_B^-1
Eigen::MatrixXd generalInverseOUCovarMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi)
{
	int md=inputs.size();
	
	//characteristic elements of SIGMA_B_INV	
	double ri = exp(-beta);
	double ei = ri / (1.0 - ri * ri);	//off-diagonal elements
	double d1= 1.0 + ri * ei;			//corner elements
	double di = 1.0 + 2.0 * ri * ei;	//general diagonal elements
	
	Eigen::MatrixXd SIGMA_B = makeSigmaBMatrix(inputs, sigma_b2, beta, kappa, pi); //simple way to remember the original elements
	
	Eigen::MatrixXd SIGMA_B_INV=Eigen::MatrixXd::Zero(md, md);
	
	//fill off-diagonal elements by inverting the 2x2 submatrix
	for(int i=0; i<md; i++){
		if(i>0){
			//calculate the off-diagonals with a normal submatrix
			double a00=SIGMA_B(i-1, i-1);
			double a01=SIGMA_B(i-1, i);
			double a11=SIGMA_B(i,i);
				
			double ainv01=a01/(a01*a01-a00*a11);
			SIGMA_B_INV(i,i-1)=ainv01;
			SIGMA_B_INV(i-1,i)=ainv01;
		}
	}		
	
	for(int i=0; i<md; i++){
		//now calculate the diagonal elements
		double sumprod=0.0;
		if(i>0){
			sumprod+=SIGMA_B(i,i-1)*SIGMA_B_INV(i,i-1);
		}
		if(i<md-1){
			sumprod+=SIGMA_B(i,i+1)*SIGMA_B_INV(i,i+1);
		}
		SIGMA_B_INV(i,i)=(1.0-sumprod)/SIGMA_B(i,i);
	}
	
	return SIGMA_B_INV;
}

//---------------------------------------------------------------------------------------------------

//PREPARES (SIGMA_E^-1 + SIGMA_B^-1)^-1
Eigen::MatrixXd makeVarBRealizationMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e)
{
	int md=inputs.size();
	Eigen::MatrixXd I=Eigen::MatrixXd::Identity(md, md);
	Eigen::MatrixXd SIGMA_E_INV = makeSigmaEInverse(inputs, sigma_e2, kappa_e); //1.0/sigma_e2 * I;
	
	//analytic solution
	Eigen::MatrixXd SIGMA_B_INV = 	generalInverseOUCovarMatrix(inputs, sigma_b2, beta, kappa, pi);
	Eigen::MatrixXd SIGMA_EINV_plus_BINV = SIGMA_E_INV + SIGMA_B_INV;
	Eigen::MatrixXd SIGMA_EINV_plus_BINV_INV = generalInvertTridiagonal(SIGMA_EINV_plus_BINV);
	
	return SIGMA_EINV_plus_BINV_INV;
}

//---------------------------------------------------------------------------------------------------

//PREPARES (SIGMA_E + SIGMA_B)^-1
Eigen::MatrixXd makeCovarMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e, double * det)
{
	//this time not optimized at all
	int md=inputs.size();
	
	Eigen::MatrixXd SIGMA_INV;
	
	//Analytic solving
	Eigen::MatrixXd SIGMA_E_INV = makeSigmaEInverse(inputs, sigma_e2, kappa_e); //1.0 / sigma_e2 * Eigen::MatrixXd::Identity(md, md);
	Eigen::MatrixXd SIGMA_EINV_plus_BINV_INV = makeVarBRealizationMatrix(inputs, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e);
	//multiply SIGMA_EINV_plus_BINV_INV twice with SIGMA_E_INV (literally )
	for(int r=0; r<md; r++){
		double einv = 1.0 / varianceOfE(inputs[r], sigma_e2, kappa_e);
		for(int c=0; c<md; c++){
			SIGMA_EINV_plus_BINV_INV(r, c) *= einv;	//once pre-multiplication (row-wise)
			SIGMA_EINV_plus_BINV_INV(c, r) *= einv;	//then post-multiplication (col-wise)
		}
	}
	SIGMA_INV = SIGMA_E_INV - SIGMA_EINV_plus_BINV_INV;
	
	if(det){
		*det=log(SIGMA_INV.determinant());
	}
		
	return SIGMA_INV;	
}

//---------------------------------------------------------------------------------------------------

Eigen::MatrixXd inflatedVarBRealization(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e, int md)
{
	Eigen::MatrixXd result;
	int dim = inputs.size();

	if(md<dim){
		
		if(md%2==0){
			std::cout<<"[Warning]: Matrix inflation requires odd size."<<std::endl;
		}
		
		result = Eigen::MatrixXd::Zero(dim, dim);
		
		//fill it
		Eigen::MatrixXd kernel;
		std::vector<double> inps;
		std::vector<double>::iterator it;
		//UR corner
		it=inputs.begin();
		inps=std::vector<double> (it, it+md);
		kernel=makeVarBRealizationMatrix(inps, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e);
		for(int r=0; r<md; r++){
			for(int c=0; c<md; c++){
				result(r,c) = kernel(r,c);
			}
		}
		
		//LR corner
		it=inputs.end()-md;
		inps=std::vector<double> (it, it+md);
		kernel=makeVarBRealizationMatrix(inps, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e);
		for(int r=0; r<md; r++){
			for(int c=0; c<md; c++){
				result(dim-md+r,dim-md+c) = kernel(r,c);
			}
		}
		
		//diagonal
		for(int k=1; k<dim-md; k++){
			it=inputs.begin()+k;
			inps=std::vector<double> (it, it+md);
			kernel=makeVarBRealizationMatrix(inps, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e);
			for(int c=0; c<md; c++){
				result(md-c-1+k,c+k)=kernel(md-c-1,c);
			}
		}
	}
	else{
		result=makeVarBRealizationMatrix(inputs, sigma_b2, beta, kappa, pi, sigma_e2, kappa_e);
	}
	return result;
}

//---------------------------------------------------------------------------------------------------

//MAKES AN ORNSTEIN-UHLENBECK STEP WITH THE GIVEN AUTOCORRELATION AND VARIANCE
double makeOUStep(double act_val, double jump_var, double beta)
{
	//make a conditional step from act_val with sigma_b2 and beta
	return act_val * exp(-beta) + invnormdist(0.0, sqrt(jump_var));
}

//---------------------------------------------------------------------------------------------------

//MAKES AN I.I.D. NOISE STEP WITH THE GIVEN VARIANCE
double makeNoiseStep(double sigma_e2, double input, double kappa_e)
{
	return invnormdist(0.0, sqrt(varianceOfE(input, sigma_e2, kappa_e)));
}

//---------------------------------------------------------------------------------------------------

#pragma mark Plain multivariate normal utilities

//MAKES A PLAIN COVARIANCE MATRIX BASED ON DATA VECTORS
Eigen::MatrixXd covarMatrix(const std::vector< std::vector<double> > & data, int startrow)
{
	int md = data.size();
	
	Eigen::MatrixXd SIGMA = Eigen::MatrixXd::Zero(md, md);
	
	//now fill up with covariances
	for(int d1=0; d1<md; d1++){
		for(int d2=d1; d2<md; d2++){
			double covar = covariance(&(data[d1]),&(data[d2]), startrow);
			SIGMA(d1, d2) = covar;
			SIGMA(d2, d1) = covar;
		}
	}
	
	return SIGMA;
}

//---------------------------------------------------------------------------------------------------

//MAKES A PLAIN COVARIANCE MATRIX BASED ON DATA VECTORS
Eigen::MatrixXd covarMatrix2(const std::vector< std::vector<double> > & data, double maxr2)
{
	double maxR2 = maxr2<=1.0 ? maxr2 : 1.0;
	int md = data.size();
	
	Eigen::MatrixXd SIGMA = Eigen::MatrixXd::Zero(md, md);
	
	//now fill up with covariances
	for(int d1=0; d1<md; d1++){
		for(int d2=d1; d2<md; d2++){
			double covar = covariance(&(data[d1]),&(data[d2]), 0);
			double correl = correlation(&(data[d1]),&(data[d2]));
			double r2 = correl * correl;
			double sigma1 = sqrt(variance(&(data[d1])));
			double sigma2 = sqrt(variance(&(data[d2])));
			if(r2 > maxR2){
				//get rid of too high correlation
				double sign_r = correl<0.0?-1.0:1.0;
				covar = sign_r * sqrt(maxR2) * sigma1 * sigma2;
			}
			SIGMA(d1, d2) = covar;
			SIGMA(d2, d1) = covar;
		}
	}
	
	return SIGMA;
}

//---------------------------------------------------------------------------------------------------

Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd & sigma)
{
	//decompose AiCiI
	Eigen::LLT<Eigen::MatrixXd> choleskyDecompositionOfSigma (sigma);
	Eigen::MatrixXd L=choleskyDecompositionOfSigma.matrixL();
	
	return L;
}

//---------------------------------------------------------------------------------------------------

std::vector<double> multivariateNormal(const Eigen::MatrixXd & L, const std::vector<double> & stddraws, const std::vector<double> & mus)
{
	//bounds check
	int md = L.rows();
	int vd = stddraws.size();
	int mud = mus.size();
	
	std::vector<double> result;
	
	if(md!=vd || md!=mud){
		return result;
	}
	
	//make Eigen-style vector
	Eigen::VectorXd stdvec = Eigen::VectorXd::Zero(md);
	
	for(int i=0; i<md; i++){
		stdvec(i)=stddraws[i];
	}
	
	//make the multiplication
	Eigen::VectorXd res = L * stdvec;
	
	result.assign(md,0.0);
	for(int i=0; i<md; i++){
		result[i]=mus[i]+res(i);
	}
	
	return result;
}

//---------------------------------------------------------------------------------------------------

bool isfinite(const Eigen::MatrixXd & x){
	const double * dd = x.data();
	for(size_t i = 0, size = x.size(); i < size; i++){
      	double temporary = (*(dd + i));
    	if(std::isnan(temporary) || std::isinf(temporary)){
        	std::cout << "Matrix::isfinite is false at item #"<<i<<" for "<<temporary<<std::endl;
        	return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------------------------------
