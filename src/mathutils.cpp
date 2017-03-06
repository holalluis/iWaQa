/*
 *  mathutils.cpp
 *  Mathematical utility functions 
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/UTILS
 *
 */ 
 
#include "mathutils.h"

#include <functional>
#include <numeric>
#include <algorithm>
#include <stdlib.h>
#include <float.h>

#include "unixtools.h"

//========================================================================================================

#pragma mark Stat functions & Co.

// Special case LOESS smoothing for uniformly spaced data (multiple samples of y are allowed at each location)

iWQVector loess(const std::vector<iWQVector> * y, double factor)
{
	iWQVector result;
	if(!y){
		return result;
	}
	int numvectors=y->size();
	if(numvectors==0){
		return result;
	}
	int npoints = 0;
	for(int i=0; i<numvectors; i++){
		int act_size=y->at(i).size();
		if(act_size>npoints){
			npoints=act_size;
		}
	}
	if(npoints==0){
		return result;
	}
	int windowhalfwidth=(int)(factor*(double)npoints/2.0);
	if(windowhalfwidth<=0){
		return result;
	}
	if(windowhalfwidth>npoints){
		windowhalfwidth=npoints/2;
	}
	
	for(int i=0; i<npoints; i++){
		
		//find the windowhalfwidth closest points to i
		int startj=(i-windowhalfwidth>=0)?(i-windowhalfwidth):0;
		int endj=(i+windowhalfwidth<npoints)?(i+windowhalfwidth):(npoints-1);
		
		double sumweight=0.0;
		double sumweightx=0.0;
		double sumweightx2=0.0;
		double sumweighty=0.0;
		double sumweightxy=0.0;
		
		for(int j=startj; j<=endj; j++){
			double weight=pow(1-pow(((i>j)?(i-j):(j-i))/((double)windowhalfwidth),3),3);
			for(int k=0; k<numvectors; k++){
				if(j<y->at(k).size()){
					sumweight+=weight;
					sumweightx+=j * weight;
					sumweightx2+=j * j * weight;
					sumweighty+=y->at(k).at(j) * weight;
					sumweightxy+=j * (y->at(k).at(j)) * weight;
				}
			}
		}
		double denom=sumweight*sumweightx2-(sumweightx*sumweightx);
		
		double a=(sumweight*sumweightxy-sumweightx*sumweighty)/denom;
		double b=(sumweightx2*sumweighty-sumweightx*sumweightxy)/denom;
		
		result.push_back(a*i+b);
	}

	return result;
}

//========================================================================================================

// Convenience wrapper to call the above on a single sample

iWQVector loess(const iWQVector * y, double factor)
{
	std::vector<iWQVector> ym;
	ym.push_back(*y);
	return loess(&ym,factor);
}

//========================================================================================================

// Fully customizable loess routine (1 dataset Y=f(X), output domain=xDomain)
iWQVector loess(const iWQVector * X, const iWQVector * Y, const iWQVector * xDomain, double factor)
{
	int i, iMin, iMax, iPoint, iMx;
	double maxDist, SumWts, SumWtX, SumWtX2, SumWtY, SumWtXY;
	double Denom, WLRSlope, WLRIntercept, xNow;
	iWQVector yLoess;
	iWQVector weight;
	iWQVector distance;
	
	int nPts=(int)(factor * X->size());
  
	//Dim mx As Variant
  
	yLoess.assign(xDomain->size(),0);
	distance.assign(X->size(),0);
	weight.assign(X->size(),0);
	
	for(int iPoint=0; iPoint<xDomain->size(); iPoint++){
		iMin = 0;
		iMax = X->size()-1;

		xNow = xDomain->at(iPoint);

		for(i=iMin; i<=iMax; i++){
			// populate x, y, distance
			distance[i] = fabs(X->at(i)-xNow);
		}
    
		for(;;){
			// find the nPts points closest to xNow
			if(iMax - iMin <= nPts){
				break;
			}
			if(distance[iMin] > distance[iMax]){
				// remove first point
				iMin++;
			}
			else if(distance[iMin] < distance[iMax]){
				// remove last point
				iMax--;
			}
			else{
				// remove both points?
				iMin++;
				iMax--;
			}
		}
		
		// Find max distance
		maxDist = -1.0;
		for(i = iMin; i<=iMax; i++){
			if(distance[i] > maxDist){ maxDist = distance[i]; }
		}

		// calculate weights using scaled distances
		for(i = iMin; i<=iMax; i++){
			weight[i] = pow(1 - pow(distance[i] / maxDist , 3), 3);
		}

		// do the sums of squares
		SumWts = 0.0;
		SumWtX = 0.0;
		SumWtX2 = 0.0;
		SumWtY = 0.0;
		SumWtXY = 0.0;
    
		for(i = iMin; i<=iMax; i++){
			SumWts += weight[i];
			SumWtX += X->at(i) * weight[i];
			SumWtX2 += X->at(i) * X->at(i) * weight[i];
			SumWtY += Y->at(i) * weight[i];
			SumWtXY += X->at(i) * Y->at(i) * weight[i];
		}
		Denom = SumWts * SumWtX2 - SumWtX * SumWtX;

		// calculate the regression coefficients, and finally the loess value
		WLRSlope = (SumWts * SumWtXY - SumWtX * SumWtY) / Denom;
		WLRIntercept = (SumWtX2 * SumWtY - SumWtX * SumWtXY) / Denom;
		yLoess[iPoint] = WLRSlope * xNow + WLRIntercept;
	}
	return yLoess;
}

//========================================================================================================

iWQVector density(iWQVector x, iWQVector q)
{
	//Simple kernel density estimation
	// x: sample, q: evaluation positions
	// returns: the estimated values of the density function of x at q
	int nx = x.size();
	int nq = q.size();
	
	iWQVector ys, ps, yps;
	
	if(nx==0 || nq==0){
		return ys;
	}
	
	double inx = 1.0 / nx;
	double inq = 1.0 / nq;
	
	std::sort(x.begin(), x.end());
	double q75 = quantile(x, 0.75, 7, true);
	double q25 = quantile(x, 0.25, 7, true);
	
	double IQR = q75 - q25;
	double h = 1.06 * IQR / 1.34 * pow((double)nx, -0.2);
	
	if(h==0.0){
		return ys;
	}
	
	//do it for the entire dataset with 512 point coverage
	int np=512;
	double xmin=min(&x);
	double xmax=max(&x);
	double dp=(xmax-xmin)/(double)(np-1);
	ps.assign(np, 0);
	ps[0] = xmin;
	for(int i=1; i<np; i++){
		ps[i] = ps[i-1] + dp;
	}
	
	yps.assign(np, 0);
	for(int i=0; i<np; i++){
		for(int j=0; j<nx; j++){
			yps[i] += inx * dnorm((ps[i] - x[j])/h);
		}
	}
	
	//normalize area
	double integral = sum(&yps) * dp;
	if(integral==0.0){
		return ys;	//no integral
	}
	
	//printf("Integral: %lf\n", integral);
	
	for(int i=0; i<np; i++){
		yps[i] /= integral;
	}
	
	/*FILE * fdens = fopen("debug_density_p.txt","w");
	for(int i=0; i<np; i++){
		fprintf(fdens, "%lf\t%lf\n", ps[i], yps[i]);
	}
	fclose(fdens);*/
	
	//get the value at qs
	ys.assign(nq, 0);
	for(int i=0; i<nq; i++){
		double q_act = q[i];
		int x1index = (int)((q_act - xmin)/dp);
		//printf("q=%lf, x1i=%d, xmin=%;f, dp=%lf\n",q_act,x1index,xmin,dp);
		if(x1index<0){
			//extrapolation below xmin
			ys[i] = interpolate(q_act, ps[0], yps[0], ps[1], yps[1]);
		}
		else if(x1index<(np-1)){
			//interpolation
			ys[i] = interpolate(q_act, ps[x1index], yps[x1index], ps[x1index+1], yps[x1index+1]);
		}
		else{
			//extrapolation above xmax
			ys[i] = interpolate(q_act, ps[np-2], yps[np-2], ps[np-1], yps[np-1]);
		}
	}
		
	return ys;
}

//--------------------------------------------------------------------------------------------------------

double interpolate(double x, double x1, double y1, double x2, double y2)
{
	if(x1==x2 && x!=x1){
		return 0.0;
	}
	double a = (y2 - y1) / (x2 - x1);
	double b = y1 - a * x1;
	
	return a * x + b;
}

//========================================================================================================

double average(const iWQVector * x)
{
	return average(x, 0);
}

//--------------------------------------------------------------------------------------------------------
double average(const iWQVector * x, int startrow)
{
	//return std::accumulate(x->begin(), x->end(), 0.0) / (double)x->size();
	double sumx=0.0;
	int n=x->size();
	int numfaulty=startrow;
	
	for(int i=startrow; i<n; i++){
		double val=x->at(i);
		if(isnan(val)){
			numfaulty++;
		}	
		else{
			sumx+=val / (double)n;	//division moved here to prevent overflows
		}
	}
	return sumx * n / (double)(n-numfaulty); // / (double) n;
}

//========================================================================================================

//Sum
double sum(const iWQVector * x)
{
	double sumx=0.0;
	int n=x->size();
	
	for(int i=0; i<n; i++){
		double val=x->at(i);
		if(!isnan(val)){
			sumx+=val;	//division moved here to prevent overflows
		}
	}
	return sumx;
}

//========================================================================================================

//Min
double min(const iWQVector * x)
{
	double min=DBL_MAX;
	int n=x->size();
		
	for(int i=0; i<n; i++){
		double val=x->at(i);
		if(!isnan(val) && !isinf(val)){
			if(val<min){
				min=val;
			}
		}	
	}
	return min;
}

//========================================================================================================

//Max
double max(const iWQVector * x)
{
	double max=-DBL_MAX;
	int n=x->size();
		
	for(int i=0; i<n; i++){
		double val=x->at(i);
		if(!isnan(val) && !isinf(val)){
			if(val>max){
				max=val;
			}
		}	
	}
	return max;
}

//========================================================================================================

//Sum of squares
double sumsquares(const iWQVector * x)
{
	double sumsqx=0.0;
	int n=x->size();
	
	for(int i=0; i<n; i++){
		double val=x->at(i);
		if(!isnan(val)){
			sumsqx+=(val * val);
		}
	}
	return sumsqx;
}

//========================================================================================================

//Variance
double variance(const iWQVector * x)
{
	double avg=average(x);
	
	double sumsqx=0.0;
	double val;
	int n=x->size();
	int numfaulty=0;
	
	for(int i=0; i<n; i++){
		val=x->at(i);
		if(isnan(val)){
			numfaulty++;
		}
		else{
			sumsqx+=(val-avg)*(val-avg)/(n-1.0);	//division moved here to prevent overflows
		}
	}
	if(n-numfaulty<=0){
		return 0.0;
	}
	else{
		return sumsqx * n / (double)(n-numfaulty);
	}
	/*int n=x->size();
	iWQVector diff(n);
	std::transform(x->begin(), x->end(), diff.begin(), std::bind2nd(std::minus<double>(), avg));
	return std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / (double)(n - 1.0);*/
}

//========================================================================================================

//Correlation
double correlation(const iWQVector * x, const iWQVector * y)
{
	double xavg=average(x);
	double yavg=average(y);
	
	double sumdevx2=0.0;
	double sumdevy2=0.0; 
	double sumdevxdevy=0.0;

	double xval, yval, devx, devy;
	
	int nx=x->size();
	int ny=y->size();
	
	int n=(nx<ny)?nx:ny;
	
	int nfaulty=0;
	
	for(int i=0; i<n; i++){
		xval=x->at(i);
		yval=y->at(i);
		
		if(isnan(xval) || isnan(yval)){	//some of them is +/-INF
			nfaulty++;
		}
		else{
			devx=xval-xavg;
			devy=yval-yavg;
		
			sumdevx2+=devx*devx;
			sumdevy2+=devy*devy;
			sumdevxdevy+=devx*devy;
		}
	}
	
	if(sumdevx2!=0.0 && sumdevy2!=0.0){
		return sumdevxdevy / sqrt(sumdevx2 * sumdevy2);
	}
	else{
		return 0.0;
	}
}

//========================================================================================================

//Covariance
double covariance(const iWQVector * x, const iWQVector * y)
{
	return covariance(x, y, 0);
}

//-------------------------------------------------------------------------------------------------------

double covariance(const iWQVector * x, const iWQVector * y, int startrow)
{
	double xavg=average(x, startrow);
	double yavg=average(y, startrow);
	
	double sumdevxdevy=0.0;

	double xval, yval, devx, devy;
	
	int nx=x->size();
	int ny=y->size();
	
	int n=(nx<ny)?nx:ny;
	
	if(n==0){
		return 0.0;
	}
	
	int nfaulty=startrow;
	
	for(int i=startrow; i<n; i++){
		xval=x->at(i);
		yval=y->at(i);
		
		if(isnan(xval) || isnan(yval)){	//some of them is +/-INF
			nfaulty++;
		}
		else{
			devx=xval-xavg;
			devy=yval-yavg;
		
			sumdevxdevy+=devx*devy;
		}
	}
	
	if(n-nfaulty<=0){
		return 0.0;
	}
	else{
		return sumdevxdevy / (double)(n-nfaulty);
	}
}

//========================================================================================================

//Quantile
double quantile(iWQVector x, double q,  unsigned short int qtype, bool sorted)
{
	/*Adapted from the perl code of Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.*/
    
	if(!x.size()){	//return 0 for empty vectors
		return 0.0;
	}
	
	//sort the data
	if(!sorted){
		std::sort(x.begin(), x.end());
	}
	
	if(qtype<1 || qtype>9){
		//error in qtype
		return 0.0;
	}
	
    // Parameters for the Hyndman and Fan algorithm
    double abcd [9][4] = {	{0.0, 0.0, 1.0, 0.0}, // inverse empirical distrib.function., R type 1
							{0.5, 0.0, 1.0, 0.0}, // similar to type 1, averaged, R type 2
							{0.5, 0.0, 0.0, 0.0}, // nearest order statistic,(SAS) R type 3
							{0.0, 0.0, 0.0, 1.0}, // California linear interpolation, R type 4
							{0.5, 0.0, 0.0, 1.0}, // hydrologists method, R type 5
							{0.0,   1.0, 0.0, 1.0}, // mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
							{1.0,  -1.0, 0.0, 1.0}, // mode-based method,(S, S-Plus), R type 7
							{1.0/3.0, 1.0/3.0, 0.0, 1.0}, // median-unbiased ,  R type 8
							{3.0/8.0, 0.25, 0.0, 1.0}   // normal-unbiased, R type 9.
						};
 
	double a = abcd[qtype-1][0];
	double b = abcd[qtype-1][1];
	double c = abcd[qtype-1][2];
	double d = abcd[qtype-1][3];
	
	int n = x.size();
    
	double g, j;
	g = modf( a + (n + b) * q - 1.0, &j);
    if(j < 0){
        return x.at(0);
	}
	else if(j >= n){
        return x.at(n-1);
	}
	
	int jj=(int)floor(j);
    
	if(g == 0){
       return x.at(jj);
    }
	else{
       return x.at(jj) + (x.at(jj+1) - x.at(jj)) * (c + d * g);
	}
}

//========================================================================================================

void sampleConfidenceLimits(iWQVector * data, double p, double * low, double * high)
{
	if(!low || !high){
		return;
	}
	
	if(p<=0.0 || p>1.0){
		*low=DBL_MAX;
		*high=-DBL_MAX;
		return;
	}
	
	//sort the data
	std::sort(data->begin(), data->end());
	
	int ndata=data->size();
	int dof=(int)((1.0-p)*(double)ndata);	//degree of freedom to shift the interval, if 0 then MIN and MAX are returned 
	
	//now we have dof varieties of possible CI span
	double minspan=DBL_MAX;	
	int minindex=-1;
	for(int i=0; i<=dof; i++){		//record spans of different CIs with the given p
		double span=data->at(ndata-1-dof+i)-data->at(i);
		if(span<minspan){
			minindex=i;
			minspan=span;
		}
	}
	
	//report the smallest CI for p
	*low=data->at(minindex);
	*high=data->at(ndata-1-dof+minindex);
	
	return;
}

//========================================================================================================

/*
Soft maximum function with overflow-safe implementation
Modified from http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
The method relies on a special quasi-arithmetic mean or generalised f-mean
but also can be interpreted as a sigmoid smoothing
*/

double SoftMaximum(double x, double y, double k)
{
	if(k<1.0){
		k=1.0;
	}
	double maximum = k * (x>=y?x:y);
	double minimum = k * (x<=y?x:y);
	return (log1p( exp(minimum - maximum) ) + maximum) / k;
}

//========================================================================================================

//Soft threshold function
double SoftThreshold(double x, double threshold, double k)
{
	//Sigmoid "step" function at threshold with k shape factor
	double exppart = -k * (x-threshold);
	if(exppart >= -50.0 && exppart <50.0){
		return 1.0 / (1.0 + exp(exppart));
	}
	else if(exppart<-50.0){		//direct solutions for extreme values
		return 0.0;
	}
	else{
		return 1.0;
	}
}

//========================================================================================================

//renamed versions of the min/max selection
double constrain_minmax(double x, double min, double max)	//equals to min(max,max(x,min))
{
	if(x<min){
		return min;
	}
	if(x>max){
		return max;
	}
	return x;
}

double constrain_min(double x, double min)	//equals to max(x,min)
{
	if(x<min){
		return min;
	}
	return x;
}

double constrain_max(double x, double max)	//equals to min(x,max)
{
	if(x>max){
		return max;
	}
	return x;
}

//========================================================================================================

double urand()		//valid only for a single thread.
{
	return rand()/((double)(RAND_MAX+1.0));
}

//========================================================================================================

double invnormdist(double mean, double sdev)	//valid only for a single thread.
{
	double U1=urand();
	double U2=urand();
	while(U1==0.0){
		U1=urand();
	}
		
	double x=sqrt(-2*log(U1))*cos(2*M_PI*U2);
	
	return sdev*x+mean;
}

//========================================================================================================

void samplerKernel(int n, double * act_values, double * new_values, double * sdevs)	//valid only for a single thread.
{
	for(int i=0; i<n; i++){
		new_values[i]=invnormdist(act_values[i],sdevs[i]);
	}
}

//========================================================================================================

#pragma mark Distributions / random number generators

// generic class functionality to spawn random generators
iWQRandomGenerator::iWQRandomGenerator(int factor)
{
	//seed random number generator (uniquely among many threads (factor = thread id))
	struct timeval tv;
	gettimeofday(&tv,NULL);
	mSeed=tv.tv_sec*tv.tv_usec*(factor+1);
	srand(mSeed);
}

//--------------------------------------------------------------------------------------------------------

//uniform random number between [0.0, 1.0[
double iWQRandomGenerator::uniform_random()
{
	return rand_r(&mSeed)/((double)(RAND_MAX+1.0));
}

//========================================================================================================

#pragma mark Uniform distribution

iWQRandomUniformGenerator::iWQRandomUniformGenerator(double min, double max, int threadid) : iWQRandomGenerator(threadid)
{
	setMin(min);
	setMax(max);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomUniformGenerator::generate()
{
	double p = uniform_random();		
	return mMin + p * (mMax - mMin);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomUniformGenerator::logLikeli(double x)
{ 
	if(x>=mMin && x<mMax){
		return -log(fabs(mMax - mMin));
	}
	else{
		return -DBL_MAX; 
	}
} 

//--------------------------------------------------------------------------------------------------------

void iWQRandomUniformGenerator::initialize(iWQDistributionSettings settings)
{
	double min=settings["min"];
	double max=settings["max"];
	setMin(min);
	setMax(max);
}

//========================================================================================================

#pragma mark Exponential distribution

iWQRandomExpGenerator::iWQRandomExpGenerator(double mean, int threadid) : iWQRandomGenerator(threadid)
{
	mLambda = 1.0;	//to provide a safe default
	setMean(mean);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomExpGenerator::generate()
{
	double p = uniform_random();
	return -log(1.0 - p) / mLambda;
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomExpGenerator::logLikeli(double x)
{
	if(x<=0.0){
		return -DBL_MAX; 
	}
	return log(mLambda) - mLambda * x;	//logarithm of mLambda * exp(-mLambda*x);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomExpGenerator::setMean(double val)
{ 
	if(val>0.0){
		mLambda = 1.0 / val; 	
	}
	//else ignore invalid val
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomExpGenerator::initialize(iWQDistributionSettings settings)
{
	double mean=settings["mean"];
	setMean(mean);
}

//========================================================================================================

#pragma mark Normal distribution

void iWQRandomNormalGenerator::generate2numbers()
{
	double U, V;
	do{
		U=uniform_random();
		V=uniform_random();
	}while(U<=0.0 || V<=0.0 || U>=1.0 || V>=1.0);
	
	double c1 = sqrt(-2.0 * log(U));
	double c2 = 2.0 * M_PI * V;
	mR1=c1*cos(c2);
	mR2=c1*sin(c2);
	mExportedCount=0;
}

//--------------------------------------------------------------------------------------------------------

iWQRandomNormalGenerator::iWQRandomNormalGenerator(double mean, double stdev, int threadid) : iWQRandomGenerator(threadid)
{
	setMean(mean);
	setStdev(stdev);
	mExportedCount=2;
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomNormalGenerator::generate()
{
	if(mExportedCount>=2){
		generate2numbers();
	}
	mExportedCount++;
	switch(mExportedCount){
		case 1: return mR1*mStdev+mAvg;
		case 2: return mR2*mStdev+mAvg;
	}
	return 0.0; //fake value, should never reach here
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomNormalGenerator::setMean(double val)
{ 
	mAvg=val; 
}
	
//--------------------------------------------------------------------------------------------------------

void iWQRandomNormalGenerator::setStdev(double val)
{ 
	mStdev=val; 
	//pre-calculations for log likelihood 
	mLogFirstPart=-0.5*log(2.0*M_PI*val*val);
	mInverse2SigmaSquare=1.0/(2.0*val*val);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomNormalGenerator::logLikeli(double x)
{
	return mLogFirstPart-(x-mAvg)*(x-mAvg)*mInverse2SigmaSquare;
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomNormalGenerator::initialize(iWQDistributionSettings settings)
{
	double mean=settings["mean"];
	double stdev=settings["stdev"];
	
	setMean(mean);
	setStdev(stdev);
}

//========================================================================================================

#pragma mark Lognormal distribution

iWQRandomLogNormalGenerator::iWQRandomLogNormalGenerator(double mean, double stdev, int threadid) : iWQRandomNormalGenerator(threadid)
{ 
	setMean(mean);
	setStdev(stdev);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomLogNormalGenerator::generate()
{
	double val=iWQRandomNormalGenerator::generate();
	return exp(val);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomLogNormalGenerator::setMean(double val)
{
	mMu=val;
	setDistParams();
} 

//--------------------------------------------------------------------------------------------------------

void iWQRandomLogNormalGenerator::setStdev(double val)
{
	mSigma=val;
	setDistParams();
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomLogNormalGenerator::setDistParams()
{
	double sdevln=sqrt(log((mSigma/mMu)*(mSigma/mMu)+1.0));
	iWQRandomNormalGenerator::setStdev(sdevln);
	double muln=log(mMu)-0.5*sdevln*sdevln;
	iWQRandomNormalGenerator::setMean(muln);
	//printf("Setdistparams: (%lf, %lf) -> (%lf, %lf)\n",mMu,mSigma,mAvg,mStdev);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomLogNormalGenerator::logLikeli(double x)
{
	if(x<=0){
		return -DBL_MAX;
	}
	return iWQRandomNormalGenerator::logLikeli(log(x));
}

//--------------------------------------------------------------------------------------------------------
	
void iWQRandomLogNormalGenerator::initialize(iWQDistributionSettings settings)
{
	double mean=settings["mean"];
	double stdev=settings["stdev"];
	
	setMean(mean);
	setStdev(stdev);
}

//========================================================================================================

#pragma mark NumPy & other utility functions for complicated distributions

//  gamma.cpp -- computation of gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
// Returns gamma function of argument 'x'.
//
// NOTE: Returns 1e308 if argument is a negative integer or 0,
//      or if argument exceeds 171.
//

double gammax(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] = {
        1.0,
        0.5772156649015329,
       -0.6558780715202538,
       -0.420026350340952e-1,
        0.1665386113822915,
       -0.421977345555443e-1,
       -0.9621971527877e-2,
        0.7218943246663e-2,
       -0.11651675918591e-2,
       -0.2152416741149e-3,
        0.1280502823882e-3,
       -0.201348547807e-4,
       -0.12504934821e-5,
        0.1133027232e-5,
       -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
       -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
       -0.36968e-11,
        0.51e-12,
       -0.206e-13,
       -0.54e-14,
        0.14e-14};

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++) {
               ga *= i;
            }
         }
         else
            ga = 1e308;
     }
     else {
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++) {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--) {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}

//--------------------------------------------------------------------------------------------------------

/*
	large parts from https://github.com/numpy/numpy/blob/master/numpy/random/mtrand/distributions.c
	Copyright 2005 Robert Kern (robert.kern@gmail.com)
	
	The random number generator's state was replaced by a simple unsigned int seed used by rand_r
*/

/* log-gamma function to support some of these distributions. The
* algorithm comes from SPECFUN by Shanjie Zhang and Jianming Jin and their
* book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
*/

double loggamma(double x)
{
    double x0, x2, xp, gl, gl0;
    long k, n;
    
    static double a[10] = {8.333333333333333e-02,-2.777777777777778e-03,
         7.936507936507937e-04,-5.952380952380952e-04,
         8.417508417508418e-04,-1.917526917526918e-03,
         6.410256410256410e-03,-2.955065359477124e-02,
         1.796443723688307e-01,-1.39243221690590e+00};
    x0 = x;
    n = 0;
    if ((x == 1.0) || (x == 2.0))
    {
        return 0.0;
    }
    else if (x <= 7.0)
    {
        n = (long)(7 - x);
        x0 = x + n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2*M_PI;
    gl0 = a[9];
    for (k=8; k>=0; k--)
    {
        gl0 *= x2;
        gl0 += a[k];
    }
    gl = gl0/x0 + 0.5*log(xp) + (x0-0.5)*log(x0) - x0;
    if (x <= 7.0)
    {
        for (k=1; k<=n; k++)
        {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}

//a bypass for the simple random number generator by me:
double rk_double(unsigned int * state)
{
	return rand_r(state)/((double)(RAND_MAX+1.0));
}

//a bypass for the simple random normal generator by me:
double rk_gauss(unsigned int * state)
{
	double f, x1, x2, r2;

    do {
        x1 = 2.0*rk_double(state) - 1.0;
        x2 = 2.0*rk_double(state) - 1.0;
        r2 = x1*x1 + x2*x2;
    }
    while (r2 >= 1.0 || r2 == 0.0);

    /* Box-Muller transform */
    f = sqrt(-2.0*log(r2)/r2);
	
	//TODO: it does not remember any of the 2 generated numbers (should do with x2),
	// so generation is not effective
	
    return f*x2;
}

double rk_normal(unsigned int *state, double loc, double scale)
{
    return loc + scale*rk_gauss(state);
}

double rk_standard_exponential(unsigned int *state)
{
    /* We use -log(1-U) since U is [0, 1) */
    return -log(1.0 - rk_double(state));
}

double rk_exponential(unsigned int *state, double scale)
{
    return scale * rk_standard_exponential(state);
}

double rk_uniform(unsigned int *state, double loc, double scale)
{
    return loc + scale*rk_double(state);
}

double rk_standard_gamma(unsigned int *state, double shape)
{
    double b, c;
    double U, V, X, Y;

    if (shape == 1.0)
    {
        return rk_standard_exponential(state);
    }
    else if (shape < 1.0)
    {
        for (;;)
        {
            U = rk_double(state);
            V = rk_standard_exponential(state);
            if (U <= 1.0 - shape)
            {
                X = pow(U, 1./shape);
                if (X <= V)
                {
                    return X;
                }
            }
            else
            {
                Y = -log((1-U)/shape);
                X = pow(1.0 - shape + shape*Y, 1./shape);
                if (X <= (V + Y))
                {
                    return X;
                }
            }
        }
    }
    else
    {
        b = shape - 1./3.;
        c = 1./sqrt(9*b);
        for (;;)
        {
            do
            {
                X = rk_gauss(state);
                V = 1.0 + c*X;
            } while (V <= 0.0);

            V = V*V*V;
            U = rk_double(state);
            if (U < 1.0 - 0.0331*(X*X)*(X*X)) return (b*V);
            if (log(U) < 0.5*X*X + b*(1. - V + log(V))) return (b*V);
        }
    }
}

double rk_gamma(unsigned int *state, double shape, double scale)
{	
    return scale * rk_standard_gamma(state, shape);
}

double rk_beta(unsigned int *state, double a, double b)
{
    double Ga, Gb;

    if ((a <= 1.0) && (b <= 1.0))
    {
        double U, V, X, Y;
        /* Use Jonk's algorithm */

        while (1)
        {
            U = rk_double(state);
            V = rk_double(state);
            X = pow(U, 1.0/a);
            Y = pow(V, 1.0/b);

            if ((X + Y) <= 1.0)
            {
                return X / (X + Y);
            }
        }
    }
    else
    {
        Ga = rk_standard_gamma(state, a);
        Gb = rk_standard_gamma(state, b);
        return Ga/(Ga + Gb);
    }
}

double rk_standard_t(unsigned int * state, double df)
{
    double N, G, X;

    N = rk_gauss(state);
    G = rk_standard_gamma(state, df/2);
    X = sqrt(df/2)*N/sqrt(G);
    return X;
}

double dse(double x, double beta)
{
  double omega_beta = pow(gammax(3.0 * (1.0+beta) / 2.0),0.5) / ((1.0+beta) * pow(gammax((1.0+beta)/2.0),1.5));
  double c_beta = pow( gammax(3.0 * (1.0+beta) / 2.0)/gammax((1.0+beta)/2.0), 1.0/(1.0+beta));
  return omega_beta * exp(-c_beta * pow(fabs(x),2.0/(1.0+beta)));          
}

double sgn(double x){
	return x>=0.0 ? 1.0 : -1.0;
}
//========================================================================================================

#pragma mark Students t distribution

iWQRandomtGenerator::iWQRandomtGenerator(int threadid, int dof) : iWQRandomGenerator(threadid)
{ 
	setDof(dof); 
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomtGenerator::setDof(double val)
{
	if(val>0){
		mNu=val;
		mLogGammaPart=log(gammax(0.5*(mNu+1))/(sqrt(mNu*M_PI)*gammax(0.5*mNu)));
	}
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomtGenerator::logLikeli(double x)
{
	return mLogGammaPart + log(pow(1.0 + x*x/mNu, -0.5*(mNu+1)));
}

//--------------------------------------------------------------------------------------------------------
	
void iWQRandomtGenerator::initialize(iWQDistributionSettings settings)
{
	double dof=settings["dof"];
	
	setDof(dof);
}

double iWQRandomtGenerator::generate()
{
	return rk_standard_t(&mSeed, mNu);
}

//========================================================================================================

#pragma mark Beta distribution

iWQRandomBetaGenerator::iWQRandomBetaGenerator(double alpha, double beta, int threadid) : iWQRandomGenerator(threadid)
{ 
	setAlpha(alpha); 
	setBeta(beta); 
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomBetaGenerator::logLikeli(double x)
{
	if(x>0.0 && x<1.0){
		return mLogGammaPart + (mAlpha-1) * log(x) + (mBeta-1) * log(1.0-x);
	}
	else{
		return -DBL_MAX;
	}
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomBetaGenerator::initialize(iWQDistributionSettings settings)
{
	double alpha=settings["alpha"];
	double beta=settings["beta"];
	
	setAlpha(alpha);
	setBeta(beta);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomBetaGenerator::setAlpha(double val)
{ 
	mAlpha=val; 
	updateLogGammaPart();
}

//--------------------------------------------------------------------------------------------------------
	
void iWQRandomBetaGenerator::setBeta(double val)
{ 
	mBeta=val; 
	updateLogGammaPart();
}	

//--------------------------------------------------------------------------------------------------------

void iWQRandomBetaGenerator::updateLogGammaPart()
{
	mLogGammaPart=log(gammax(mAlpha+mBeta)/(gammax(mAlpha)*gammax(mBeta)));
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomBetaGenerator::generate()
{
	return rk_beta(&mSeed, mAlpha, mBeta);
}

//========================================================================================================

#pragma mark Gamma distribution

iWQRandomGammaGenerator::iWQRandomGammaGenerator(double k, double theta, int threadid) : iWQRandomGenerator(threadid)
{ 
	setK(k); 
	setTheta(theta); 
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomGammaGenerator::logLikeli(double x)
{
	if(x>0.0){
		return (mK-1.0)*log(x) - x / mTheta + mLogGammaPart;
	}
	else{
		return -DBL_MAX;
	}
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomGammaGenerator::initialize(iWQDistributionSettings settings)
{
	double k=settings["k"];
	double theta=settings["theta"];
	
	setK(k);
	setTheta(theta);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomGammaGenerator::setK(double val)
{ 
	if(val > 0.0){
		mK=val; 
		updateLogGammaPart();
	}
}

//--------------------------------------------------------------------------------------------------------
	
void iWQRandomGammaGenerator::setTheta(double val)
{ 
	if(val > 0.0){
		mTheta=val; 
		updateLogGammaPart();
	}
}	

//--------------------------------------------------------------------------------------------------------

void iWQRandomGammaGenerator::updateLogGammaPart()
{
	mLogGammaPart= -loggamma(mK) - mK * log(mTheta);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomGammaGenerator::generate()
{
	return rk_gamma(&mSeed, mK, mTheta);
}

//========================================================================================================

#pragma mark SEP distribution

iWQRandomSEPGenerator::iWQRandomSEPGenerator(double beta, double xi, int threadid) : iWQRandomGenerator(threadid)
{ 
	setBeta(beta); 
	setXi(xi); 
	setMu(0.0);
	setSigma(1.0);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomSEPGenerator::logLikeli(double x)
{
	double stdx = (x - mMu)/mSigma;
	double newx = pow(mXi,-sgn(mMuXi + mSigmaXi * stdx)) * (mMuXi + mSigmaXi * stdx); //newx <- xi^(-sign(mu_xi + sigma_xi * x)) * (mu_xi + sigma_xi * x)
	double dse = mOmegaBeta * exp(-mCBeta * pow(fabs(newx),2.0/(1.0+mBeta)));
	return log(2.0 / mSigma * mSigmaXi / (mXi + 1.0/mXi) * dse);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomSEPGenerator::initialize(iWQDistributionSettings settings)
{
	double beta=settings["beta"];
	double xi=settings["xi"];
	double mu=settings["mean"];
	double sigma=settings["stdev"];
	
	setBeta(beta);
	setXi(xi);
	setMu(mu);
	setSigma(sigma);
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomSEPGenerator::setBeta(double val)
{ 
	mBeta = val;
	updateSkewPart(); 
}

//--------------------------------------------------------------------------------------------------------
	
void iWQRandomSEPGenerator::setXi(double val)
{ 
	if(val>0.0){
		mXi=val;
		updateSkewPart(); 
	}
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomSEPGenerator::setSigma(double val)
{ 
	if(val>0.0){
		mSigma=val;
	}
}

//--------------------------------------------------------------------------------------------------------

void iWQRandomSEPGenerator::updateSkewPart()
{
	//SE parts
	mOmegaBeta = pow(gammax(3.0 * (1.0+mBeta) / 2.0),0.5) / ((1.0+mBeta) * pow(gammax((1.0+mBeta)/2.0),1.5));
	mCBeta = pow(gammax(3.0 * (1.0+mBeta) / 2.0)/gammax((1.0+mBeta)/2.0), 1.0/(1.0+mBeta));
	//skew parts
	mM1 = gammax(1.0+mBeta) / (pow(gammax(3.0*(1.0+mBeta)/2.0),0.5)*pow(gammax((1.0+mBeta)/2.0),0.5));
	mM2 = 1.0;
	mMuXi = mM1 * (mXi - 1.0 / mXi); 
	mSigmaXi = sqrt((mM2 - mM1 * mM1)*(mXi * mXi + 1.0 / (mXi * mXi)) + 2.0 * mM1 * mM1 - mM2);
}

//--------------------------------------------------------------------------------------------------------

double iWQRandomSEPGenerator::generate()
{
	//return rk_gamma(&mSeed, mK, mTheta);
	//1. generate a standard gamma number
	double gt = rk_gamma(&mSeed, (1.0+mBeta)/2.0, 1.0); //gt <- rgamma(n, shape=(1+beta)/2, scale=1)
	//2. generate a random sign with uniform probability
	double st = -1.0 + (urand()<0.5 ? 2.0 : 0.0); //st <- -1 + 2 * as.numeric(runif(n)<0.5)
	//3. get the standard SE sample
	double EPt = st * pow(fabs(gt), (1.0 + mBeta)/2.0) * pow(gammax((1.0+mBeta)/2.0),0.5)/pow(gammax(3.0*(1.0+mBeta)/2.0),0.5); //EPt <- st * abs(gt)^((1+beta)/2) * (gamma((1+beta)/2)^(1/2))/(gamma(3*(1+beta)/2)^(1/2))
	//4. generate a random sign with p=1 - xi/(xi + 1/xi)
	double plim = mXi / (mXi + 1.0 / mXi); //plim <- xi/(xi + 1/xi)
	double wt = -1 + (urand()<plim ? 2.0 : 0.0); //wt <- -1 + 2 * as.numeric(runif(n)<plim)
	//5. compute the nonstandard SEP sample
	double SEPt = - wt * fabs(EPt) * pow(mXi,wt); //SEPt <- - wt * abs(EPt) * (xi^wt)
	//6. standardize the sample
	double at = -(SEPt + mMuXi)/mSigmaXi; //at <- -(SEPt + mu_xi)/sigma_xi
	//ready: return transformed value
	return at * mSigma + mMu;
}

//--------------------------------------------------------------------------------------------------------

#pragma mark Snippets for a future truncated normal distribution

//Snippets for a future proper truncated normal distribution

//BEGIN PART from http://www.johndcook.com/normal_cdf_inverse.html

// compute log(1+x) without losing precision for small values of x
double LogOnePlusX(double x)
{
    if (x <= -1.0)
    {
        /*std::stringstream os;
        os << "Invalid input argument (" << x 
           << "); must be greater than -1.0";
        throw std::invalid_argument( os.str() );*/
		return DBL_MAX;	//we don't want exceptions here
    }

    if (fabs(x) > 1e-4)
    {
        // x is large enough that the obvious evaluation is OK
        return log(1.0 + x);
    }

    // Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
    // Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
    return (-0.5*x + 1.0)*x;
}

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        /*std::stringstream os;
        os << "Invalid input argument (" << p 
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );*/
		return DBL_MAX;	//we don't want exceptions here
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}
//END PART

//simple standard normal CDF
double pnorm(double x)
{
	return 0.5 * (1.0 + erf(x/1.414213562));
}

//simple standard normal CDF
double lpnorm(double x)
{
	if(x < -4.0){
		return (log(0.5) - x * x / 2 + log(erfc(-x / 1.414213562)) - 2.0) / 1.964;
	}
	return log(pnorm(x));
}

//random generation from a truncated normal
double rtnorm(double mean, double sd, double * lower, double * upper)
{
	double Fi_a = 0.0;
	double Fi_b = 1.0;
	if(lower){
		double a = (*lower - mean)/sd;
		Fi_a = pnorm(a);
	}
	if(upper){
		double b = (*upper - mean)/sd;
		Fi_b = pnorm(b);
	}
	double p = Fi_a + urand() * (Fi_b - Fi_a);
	
	double xstd = NormalCDFInverse(p);
	
	double x = sd * xstd + mean;
	return x;
}

//normal density function
double dnorm(double x)
{
	return 1.0 / sqrt(2.0 * M_PI) /*M_1_SQRT_2PI*/ * exp(-0.5 * x * x);
}

//========================================================================================================

#pragma mark Nash-Sutcliffe index

double NS(std::vector<double> * measured, std::vector<double> * modelled, int start, int end)
{
	//check dimensions
	if(!measured || !modelled){
		return DBL_MAX;
	}
	if(measured->size()!=modelled->size()){
		return DBL_MAX;
	}
	int n=measured->size();
	if(n==0){
		return 0.0;
	}
	if(end==-1){
		end=n;
	}
	if(start>=n || start>=end){
		return DBL_MAX;
	}
	double averageme = 0.0;
	for(int i=start; i<end; i++){
		averageme+=measured->at(i);
	}
	averageme /= ((double)(end-start));
	double sumsqdev=0.0;
	double sumsqvar=0.0;
	for(int i=start; i<end; i++){
		sumsqdev+=pow(modelled->at(i) - measured->at(i), 2.0);
		sumsqvar+=pow(measured->at(i) - averageme, 2.0);
	}
	return sumsqdev / sumsqvar;
}

//========================================================================================================

#pragma mark Box-Cox transformation

double boxcox_transform(double lambda_1, double lambda_2, double value, bool * error)
{
	//same as in NSBoxCox
	double result=0.0;
	double biased=value+lambda_2;
	if(lambda_1==1.0){	//shortcut for no transform
		return biased;
	}
	if(biased<=0.0){
		if(error){
			*error=true;
		}
		biased=DBL_MIN;
	}
	
	if(lambda_1!=0.0){
		//general case
		result=(pow(biased,lambda_1)-1)/lambda_1;
	}
	else{
		result=log(biased);
	}
	
	return result;
}

//--------------------------------------------------------------------------------------------------------

double boxcox_retransform(double lambda_1, double lambda_2, double value, bool * error)
{
	//backward BoxCox
	double result=0.0;
	if(lambda_1==1.0){	//shortcut for no transform
		return value - lambda_2;
	}
	if(lambda_1!=0.0){
		//general case
		double base = lambda_1 * value + 1.0;
		double exponent = 1.0/lambda_1;
		if((base < 0.0 && floor(exponent)!=exponent) || (base == 0.0 && exponent < 0.0)){
			if(error){
				*error=true;
			}
			return result;
		}
		result=pow(base, exponent)-lambda_2;
	}
	else{
		result=exp(value)-lambda_2;
	}
	return result;
}

//========================================================================================================

#pragma mark Matrix allocation / deletion utilities

double ** AllocMatrix(int size)
{
	if(size<=0){
		return NULL;
	}
	double ** result=new double * [size];
	for(int i=0; i<size; i++){
		result[i]=new double [size];
		for(int j=0; j<size; j++){
			result[i][j]=0.0;
		}
	}
	return result;
}

//--------------------------------------------------------------------------------------------------------

void DeleteMatrix(double ** m, int size)
{
	if(m==NULL){
		return;
	}
	for(int i=0; i<size; i++){
		delete [] m[i];
	}
	delete [] m;
}

