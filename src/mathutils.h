/*
 *  mathutils.h
 *  Mathematical utility functions 
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/UTILS
 *
 */ 
 
#include <vector>
#include <map>
#include <string>
#include <math.h>

#ifndef mathutils_h
#define mathutils_h

// Mathematical utilities 

//	t, Beta and Gamma random number generators were adopted from the numPy module
//	https://github.com/numpy/numpy/blob/master/numpy/random/mtrand/distributions.c
//	Copyright 2005 Robert Kern (robert.kern@gmail.com)

//----------------------------------------------------------------------------------

// General vector type(def)
typedef std::vector<double> iWQVector;

//----------------------------------------------------------------------------------

#pragma mark simple statistical functions

// LOESS local regression

iWQVector loess(const iWQVector * y, double factor=0.25);
iWQVector loess(const std::vector<iWQVector> * y, double factor = 0.25);
iWQVector loess(const iWQVector * x, const iWQVector * y, const iWQVector * xout, double factor = 0.25);

// Simple kernel density estimation
iWQVector density(iWQVector x, iWQVector q);

//linear interpolators
double interpolate(double x, double x1, double y1, double x2, double y2);

// Average
double average(const iWQVector * x);
double average(const iWQVector * x, int startrow);

// Sum
double sum(const iWQVector * x);

// Variance
double variance(const iWQVector * x);

// Correlation
double correlation(const iWQVector * x, const iWQVector * y);

//Covariance
double covariance(const iWQVector * x, const iWQVector * y);
double covariance(const iWQVector * x, const iWQVector * y, int startrow);	//allow calculation for a subset too

// Quantile
double quantile(iWQVector x, double q,  unsigned short int qtype=7, bool sorted=false);

//Sum of squares
double sumsquares(const iWQVector * x);

//Min
double min(const iWQVector * x);

//Max
double max(const iWQVector * x);

//Asymmetric confidence interval estimation from sample
void sampleConfidenceLimits(iWQVector * data, double p, double * low, double * high);

//Soft maximum function
double SoftMaximum(double x, double y, double k);

//Soft threshold function
double SoftThreshold(double x, double threshold, double k);

//easy parameter constraining functions
double constrain_minmax(double x, double min, double max);	//equals to min(max,max(x,min))
double constrain_min(double x, double min);					//equals to max(x,min)
double constrain_max(double x, double max);					//equals to min(x,max)

//----------------------------------------------------------------------------------

#pragma mark Quick & dirty random number generators (for 1st phase of MCMC)

//Uniform random numbers between 0 and 1
double urand();

//Quick inverse normal generation with dynamic mean and stdev
double invnormdist(double mean, double sdev);

//Normally distributed sampler kernel for MCMC
void samplerKernel(int n, double * act_values, double * new_values, double * sdevs);

//Random number generation for a truncated normal distribution
double rtnorm(double mean, double sd, double * lower, double * upper);

//Std normal CDF
double pnorm(double x);

//Log of std normal CDF
double lpnorm(double x);

//Std normal PDF
double dnorm(double x);

//----------------------------------------------------------------------------------

#pragma mark Fully functional random number generators & distributions 

typedef std::map<std::string, double> iWQDistributionSettings; //plain data structure for unified parameter initialization

// Partially abstract base class for random number generators (distributions)
class iWQRandomGenerator
{
public:
		iWQRandomGenerator(int factor);									//this initializes the random number generator 
																		//and stores mSeed (needed for multithreading) 
		virtual ~iWQRandomGenerator(){ }
		virtual double generate()=0;									//forward operation: get a random number
		virtual double logLikeli(double x)=0;							//backward operation: get probability/likelihood
		virtual void initialize(iWQDistributionSettings settings)=0;	//uniform parameter initializer method
protected:
		unsigned int mSeed;
		double uniform_random();
};

typedef iWQRandomGenerator iWQDistribution;								//to make it more clear that they are actually distributions

//--------------------------------------------------------------------------------------------------------

//Constrained uniform distribution
class iWQRandomUniformGenerator : public iWQRandomGenerator
{
protected:
	double mMin;
	double mMax;
public:
	iWQRandomUniformGenerator(double min=0.0, double max=1.0, int threadid=0);
	virtual double generate();
	virtual double min(){ return mMin; }
	virtual void setMin(double val){ mMin=val; } 
	virtual double max(){ return mMax; }
	virtual void setMax(double val){ mMax=val; } 
	virtual double logLikeli(double x); 
	virtual void initialize(iWQDistributionSettings settings);
};

//--------------------------------------------------------------------------------------------------------

// Exponential distribution & random number generator
class iWQRandomExpGenerator : public iWQRandomGenerator
{
protected:
	double mLambda;
public:
	iWQRandomExpGenerator(double mean=1.0, int threadid=0);
	virtual double generate();
	virtual double mean(){ return 1.0 / mLambda; }
	virtual void setMean(double val);
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//--------------------------------------------------------------------------------------------------------

// Normal distribution & random number generator (Box-Mueller method)
class iWQRandomNormalGenerator : public iWQRandomGenerator
{
protected:
	double mStdev;
	double mAvg;
		
private:
	double mR1;
	double mR2;
	int mExportedCount;

	void generate2numbers();

	//for fast likelihood analysis
	double mLogFirstPart;
	double mInverse2SigmaSquare;	
	
public:
	iWQRandomNormalGenerator(double mean=0.0, double stdev=1.0, int threadid=0);
	virtual double generate();
	virtual double mean(){ return mAvg; }
	virtual double stdev(){ return mStdev; }
	virtual void setMean(double val);
	virtual void setStdev(double val);
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//--------------------------------------------------------------------------------------------------------

// Lognormal distribution & random number generator (based on the Normal case)
class iWQRandomLogNormalGenerator : public iWQRandomNormalGenerator
{
protected:
	double mMu;
	double mSigma;
private:
	void setDistParams();
public:
	iWQRandomLogNormalGenerator(double mean=1.0, double stdev=1.0, int threadid=0);
	virtual double generate();
	virtual void setMean(double val);
	virtual void setStdev(double val);
	virtual double mean(){ return mMu; }
	virtual double stdev(){ return mSigma; }
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//--------------------------------------------------------------------------------------------------------

//t distribution
class iWQRandomtGenerator : public iWQRandomGenerator
{
protected:
	double mNu;
	
private:
	//for fast likelihood analysis
	double mLogGammaPart;

public:
	iWQRandomtGenerator(int threadid=0, int dof=1);
	virtual double generate();
	virtual double dof(){ return mNu; }
	virtual void setDof(double val);
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//----------------------------------------------------------------------------------

//beta distribution
class iWQRandomBetaGenerator : public iWQRandomGenerator
{
protected:
	double mAlpha;
	double mBeta;

private:	
	//for fast likelihood analysis
	double mLogGammaPart;
	void updateLogGammaPart();	
	
public:
	iWQRandomBetaGenerator(double alpha=1.0, double beta=1.0, int threadid=0);
	virtual double generate();	
	virtual double alpha(){ return mAlpha; }
	virtual void setAlpha(double val);
	virtual double beta(){ return mBeta; }
	virtual void setBeta(double val);
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//----------------------------------------------------------------------------------

//gamma distribution
class iWQRandomGammaGenerator : public iWQRandomGenerator
{
protected:
	double mK;
	double mTheta;

private:	
	//for fast likelihood analysis
	double mLogGammaPart;
	void updateLogGammaPart();	
	
public:
	iWQRandomGammaGenerator(double k=1.0, double theta=1.0, int threadid=0);
	virtual double generate();	
	virtual double k(){ return mK; }
	virtual void setK(double val);
	virtual double theta(){ return mTheta; }
	virtual void setTheta(double val);
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//----------------------------------------------------------------------------------

//SEP distribution
class iWQRandomSEPGenerator : public iWQRandomGenerator
{
protected:
	double mBeta;
	double mXi;
	double mMu;
	double mSigma;
private:	
	double mM1;
	double mM2;
	double mMuXi;
	double mSigmaXi;
	double mCBeta;
	double mOmegaBeta;
	void updateSkewPart();
public:
	iWQRandomSEPGenerator(double beta=0.0, double xi=1.0, int threadid=0);
	virtual double generate();	
	virtual double beta(){ return mBeta; }
	virtual void setBeta(double val);
	virtual double xi(){ return mXi; }
	virtual void setXi(double val);
	virtual void setMu(double val){ mMu = val; }
	virtual double mu(){ return mMu; }
	virtual void setSigma(double val);
	virtual double sigma(){ return mSigma; }
	virtual double logLikeli(double x);
	virtual void initialize(iWQDistributionSettings settings);
};

//----------------------------------------------------------------------------------

#pragma mark Convenience generic matrix allocation

// Convenience matrix allocators (for simple storage)
double ** AllocMatrix(int size);
void DeleteMatrix(double ** m, int size);

//----------------------------------------------------------------------------------

#pragma mark Hydrology indices

// Upside-down Nash-Sutcliffe index
double NS(std::vector<double> * measured, std::vector<double> * modelled, int start=0, int end=-1);

//----------------------------------------------------------------------------------

#pragma mark Power transformations

// Box-Cox transformation
double boxcox_transform(double lambda_1, double lambda_2, double value, bool * error);
double boxcox_retransform(double lambda_1, double lambda_2, double value, bool * error);

#endif
