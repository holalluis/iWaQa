/*
 *  biasmatrices.h
 *  Helper functions for complicated likelihood calculations
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */

#include <vector>
#include "Eigen/Dense"

#ifndef biasmatrices_h
#define biasmatrices_h

//---------------------------------------------------------------------------------------------------

//RETURNS THE JUMP VARIANCE OF B
double jumpVarianceOfB(double sigma_b2, double beta, double kappa, double pi, double input);

//RETURNS THE VARIANCE OF E
double varianceOfE(double input, double sigma_e2, double kappa_e);

//PREPARES SIGMA_B
Eigen::MatrixXd makeSigmaBMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi);

//PREPARES SIGMA_E
Eigen::MatrixXd makeSigmaEMatrix(std::vector<double> inputs, double sigma_e2, double kappa_e);

//INVERTS ANY TRIDIAGONAL MATRIX
Eigen::MatrixXd generalInvertTridiagonal(Eigen::MatrixXd & T);

//SUPPLIES SIGMA_B^-1
Eigen::MatrixXd generalInverseOUCovarMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi);

//SUPPLIES SIGMA_E^-1
Eigen::MatrixXd makeSigmaEInverse(std::vector<double> inputs, double sigma_e2, double kappa_e);

//PREPARES (SIGMA_E^-1 + SIGMA_B^-1)^-1
Eigen::MatrixXd makeVarBRealizationMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e);

//PREPARES (SIGMA_E + SIGMA_B)^-1
Eigen::MatrixXd makeCovarMatrix(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e, double * det);

//RETURNS A FULL-SIZE (SIGMA_E^-1 + SIGMA_B^-1)^-1
Eigen::MatrixXd inflatedVarBRealization(std::vector<double> inputs, double sigma_b2, double beta, double kappa, double pi, double sigma_e2, double kappa_e, int md);

//MAKES AN ORNSTEIN-UHLENBECK STEP WITH THE GIVEN AUTOCORRELATION AND VARIANCE
double makeOUStep(double act_val, double jump_var, double beta);

//MAKES AN I.I.D. NOISE STEP WITH THE GIVEN VARIANCE
double makeNoiseStep(double sigma_e2, double input, double kappa_e);

//---------------------------------------------------------------------------------------------------

//MAKES A PLAIN COVARIANCE MATRIX BASED ON DATA VECTORS
Eigen::MatrixXd covarMatrix(const std::vector< std::vector<double> > & data, int startrow=0);

//MAKES A PLAIN COVARIANCE MATRIX BASED ON DATA VECTORS
Eigen::MatrixXd covarMatrix2(const std::vector< std::vector<double> > & data, double maxr2);

//PERFORMS THE CHOLESKY DECOMPOSITION ON A COVARIANCE MATRIX
Eigen::MatrixXd choleskyDecomposition(Eigen::MatrixXd & sigma);

//MAKES A MULTIVARIATE DRAW FROM INDEPENDENT STD NORMAL NUMBERS
std::vector<double> multivariateNormal(const Eigen::MatrixXd & L, const std::vector<double> & stddraws, const std::vector<double> & mus);

//---------------------------------------------------------------------------------------------------

//NA and INF checks for matrices
bool isfinite(const Eigen::MatrixXd & x);

#endif
