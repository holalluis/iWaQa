/*
 *  evaluator.cpp
 *  General likelihood manager
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/LIKELIHOOD
 *
 */
 
#include "evaluator.h"
#include "model.h"
#include "solver.h"
#include "datatable.h"
#include "evaluatormethod.h"
#include "particleswarm.h"
#include "mathutils.h"
#include "filter.h"
#include "script.h"

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#pragma mark Evaluator

//-----------------------------------------------------------------------------------

iWQEvaluator::iWQEvaluator()
{
	mSolver=0;
	mCommonParameters=0;
	mDataTable=0;
	mEvaluatorMethods.clear();
	mEvaluatorWeights.clear();
	printWarnings=true;
	returnUnstableSolutions=false;
	
	//event-based services
	mEvaluateStartRow=-1;	//not specified
	mEvaluateEndRow=-1;		//not specified
	mModelState.clear();
	mRainColPtr=NULL;
	
	//temporary PSO params
	PSOActive=false;
	PSOMaxNumRounds=100;
	PSOMaxIdleRounds=10;
	PSOSwarmSize=20;
	
	//temporary NMS params
	NMSActive=true;
	NMSMaxNumRounds=100;
	NMSTolerance=1E-7;

}

//-----------------------------------------------------------------------------------

iWQEvaluator::~iWQEvaluator()
{
	//dispose evaluator methods
	for(int i=0; i<mEvaluatorMethods.size(); i++){
		if(mEvaluatorMethods[i]){
			delete mEvaluatorMethods[i];
		}
	}
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setDataTable(iWQDataTable * datatable)
{
	if(datatable->timePort()){
		mDataTable=datatable;
	}
	else{
		printf("[Error]: Data table must have a TIME field.\n");
		mDataTable=0;
	}
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setComparisonLinks(iWQComparisonLinkSet links)
{
	mComparisonLinks=links;
}

//-----------------------------------------------------------------------------------

bool iWQEvaluator::setPredictiveMode(bool mode)
{
	for(int i=0; i<mComparisonLinks.size(); i++){
		mComparisonLinks[i].setPredictiveMode(mode);
		
		//propagate this through the evaluator methods
		for(int i=0; i<mEvaluatorMethods.size(); i++){
			if(mEvaluatorMethods[i]){
				mEvaluatorMethods[i]->setLinkPredictiveMode(mode);
			}
		}
	}
	return true;
}

//-----------------------------------------------------------------------------------
	
void iWQEvaluator::setSolver(iWQSolver * solver)
{
	mSolver=solver;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setFilters(std::vector<iWQFilter *> filters)
{
	mFilters = filters;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setPreScripts(std::vector<iWQScript> pres)
{
	mPreScripts = pres;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setPostScripts(std::vector<iWQScript> posts)
{
	mPostScripts = posts;
}
	
//-----------------------------------------------------------------------------------

void iWQEvaluator::setParameters(iWQParameterManager * parameters)
{
	mCommonParameters=parameters;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setInitialValues(iWQInitialValues * inits)
{
    mInitVals=inits;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setEvaluatorMethods(std::vector<iWQEvaluatorMethod *> methods)
{
	//empty existing methods AND their weights (if any)
	for(int i=0; i<mEvaluatorMethods.size(); i++){
		if(mEvaluatorMethods[i]){
			delete mEvaluatorMethods[i];
		}
	}
	mEvaluatorMethods.clear();
	mEvaluatorWeights.clear();
	
	//load the new ones with a default weight of 1.0
	mEvaluatorMethods=methods;
	for(int i=0; i<mEvaluatorMethods.size() && i<mComparisonLinks.size(); i++){
		if(mEvaluatorMethods[i]){
			//wire up general connections and set default weight
			mEvaluatorMethods[i]->setDataTable(mDataTable);
			mEvaluatorMethods[i]->setComparisonLink(mComparisonLinks[i]);
			mEvaluatorWeights.push_back(1.0);
		}
	}
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::setEvaluatorWeights(std::vector<double> weights)
{
	for(int i=0; i<mEvaluatorWeights.size() && i<weights.size(); i++){
		mEvaluatorWeights[i]=weights[i];
	}
}

//-----------------------------------------------------------------------------------

double iWQEvaluator::evaluate(double * values, int numpars)
{
	//direct method to speed up NelderMead
	if(mCommonParameters){
		mCommonParameters->setPlainValues(values, numpars);
		return evaluate();
	}
	else{
		printf("[Error]: Evaluator was misconfigured.\n");
		return -DBL_MAX;
	}
}

//-----------------------------------------------------------------------------------

double iWQEvaluator::evaluate(std::vector<double> values)
{
	//need to update the parameters 1st
	if(mCommonParameters){
		mCommonParameters->setPlainValues(values);
		return evaluate();
	}
	else{
		printf("[Error]: Evaluator was misconfigured.\n");
		return -DBL_MAX;
	}
}

//-----------------------------------------------------------------------------------

double iWQEvaluator::evaluate()
{
	if(!mDataTable || !mSolver || !mInitVals || !mCommonParameters || mComparisonLinks.size()==0 || mDataTable->timePort()==NULL || mEvaluatorMethods.size()==0 || mComparisonLinks.size()!=mEvaluatorMethods.size() || mEvaluatorWeights.size()!=mEvaluatorMethods.size()){
		printf("[Error]: Evaluator was misconfigured.\n");
		return DBL_MAX;	
	}
		
	//pre-filter: don't run the model when prior likelihood is invalid
	//TODO	
		
	int startrow = (mEvaluateStartRow!=-1)?mEvaluateStartRow:0;
	int endrow = (mEvaluateEndRow!=-1)?mEvaluateEndRow:mDataTable->numRows();
	
	if(mRainColPtr!=NULL){
		//populate inputs from the parameters
		std::vector<double> plainparamvals = mCommonParameters->plainValues();
		std::vector<std::string> plainparamnames = mCommonParameters->namesForPlainValues();
		int numpars = plainparamnames.size();
		for(int x=0; x<numpars; x++){
			std::stringstream s (plainparamnames[x].substr(1));
			int index;
			s >> index;
			index+=startrow;
			if(index>=0 && index<mDataTable->numRows()){
				mDataTable->setRow(index);
				*mRainColPtr = plainparamvals[x];
			}
		}
	}
	
	//run PRE scripts
	bool scriptsok=true;
	for(int s=0; s<mPreScripts.size(); s++){
		bool thisok = mPreScripts[s].execute();
		if(!thisok){
			printf("[Error]: Script \"%s\" failed to run correctly (return code=%d).\n",mPreScripts[s].commandString().c_str(),mPreScripts[s].returnStatus());
			scriptsok=false;
		}
	}
	
	//run the model
	mDataTable->setRow(startrow);	//step where the state is given
	double * t=mDataTable->timePort();
	double prev_t = *t;
	bool stable=true;
	iWQInitialValues * yfeed=NULL;
	if(startrow==0){
		yfeed=mInitVals;
	}
	else{
		if(mModelState.size()==0){
			printf("[Error]: Model state is empty for starting point of partial run.\n");
		}
		else{
			mSolver->setModelState(mModelState);
		}
	}
	
	//run the model from startrow (step right after invocation, so real output comes from startrow+1)	
	int firsterrorrow = -1; 		//the first row where instability occurs
	double firsterrort = -DBL_MAX;
	std::vector<iWQModel *> wrongs;
	wrongs.clear();
	for(int row=startrow; row<endrow; row++){
		if(!mDataTable->stepRow()){
			printf("[Error]: Partial run stepped beyond the end of data table.\n");
			break;
		}
		if(!mSolver->solve1Step(prev_t, *t, yfeed)){
			stable=false;
			firsterrorrow = mDataTable->pos();
			firsterrort = prev_t;
			wrongs = mSolver->modelsThatDidNotSolve();
		}
		prev_t = *t;
		yfeed=NULL;
	}
	
	//run is over
	
	//run POST scripts
	for(int s=0; s<mPostScripts.size(); s++){
		bool thisok = mPostScripts[s].execute();
		if(!thisok){
			printf("[Error]: Script \"%s\" failed to run correctly (return code=%d).\n",mPostScripts[s].commandString().c_str(),mPostScripts[s].returnStatus());
			scriptsok=false;
		}
	}
	
	//activate filters
	for(int i=0; i<mFilters.size(); i++){
		iWQFilter * f = mFilters[i];
		if(f){
			f->filter();
		}
	}
	
	//check the stability of the solution
	if(!stable){
		if(printWarnings){
			printf("[Warning]: Numerical stability could not be achieved with the minimal stepsize of %e.\n",mSolver->minStepLength());
			std::vector<std::string> parnames=mCommonParameters->namesForPlainValues();
			std::vector<double> parvalues=mCommonParameters->plainValues();
			for(int i=0; i<parnames.size(); i++){
				printf("%s=%g  ",parnames[i].c_str(), parvalues[i]);
			}
			printf("\n");
			if(wrongs.size()){
				sort(wrongs.begin(), wrongs.end() );
				wrongs.erase( unique( wrongs.begin(), wrongs.end() ), wrongs.end() );	//eliminate duplicates in wrongs
				printf("*** Models causing this error ***\n");
				for(int i=0; i<wrongs.size(); i++){
					printf("\t%s (%s)\n", wrongs[i]->modelId().c_str(), wrongs[i]->modelType().c_str());
				}
			}
			else{
				printf("Strange: Despite the error there are no faulty models reported.\n");
			}
			if(firsterrorrow!=-1){
				printf("*** Error location ***\n");
				printf("\trow: #%d\n",firsterrorrow);
				printf("\tstarting time coordinate: %lf\n",firsterrort);
			}
			else{
				printf("Strange: Despite the error there is no location reported.\n");
			}
		}
		if(!returnUnstableSolutions){
			return DBL_MAX;
		}
	}
	
	if(!scriptsok){
		printf("[Warning]: At least one of the scripts did not execute properly.\n");
		if(!returnUnstableSolutions){
			return DBL_MAX;
		}	
	}

	//calculate normal evaluation result
	//now collect all results from the evaluator methods
	double result=0.0;
	bool needspriors=mEvaluatorMethods[0]->priorsApply();
	if(needspriors && mCommonParameters){
		double priorlikelihood = mCommonParameters->logLikelihood(printWarnings);
		result -= priorlikelihood;		//needs to be reversed for minimisation
	}
	if(result > -DBL_MAX){ 	//no need to evaluate when prior likelihood is -INF
		for(int i=0; i<mEvaluatorMethods.size(); i++){
			mEvaluatorMethods[i]->updateDynamicParams();		//need to refresh the evaluation parameters
			double newres=mEvaluatorMethods[i]->evaluate(startrow+1,endrow);
			result+=mEvaluatorWeights[i]*newres; 	//they return negative values for minimisation
			if(printWarnings && (isnan(newres) || isinf(newres) || newres== DBL_MAX || newres== -DBL_MAX)){
				printf("[Warning]: Log likelihood of %s = %g\n",mEvaluatorMethods[i]->modelFieldName().c_str(),newres);
			}
		}
	}
	
	if(isnan(result) || isinf(result)){  //final check, so that no NaN or Inf gets outside
		result=DBL_MAX;
	}
	return result;	//need to evaluate from startrow+1 because startrow holds the initial values, which are sometimes missing
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::sequentialCalibrateParameters(std::string eventflagfield)
{
	//check if eventflagfield points to a valid column
	if(!mDataTable->portForColumn(eventflagfield)){
		printf("[Error]: There is no data column with name: %s\n",eventflagfield.c_str());
	}
	
	//need to have columns for each adjustable parameter in the data table
	std::vector<std::string> freeParams;	//data table names for the adjustable params
	std::vector<double *> paramPorts;
	
	//automatically generate free params from the common parameter storage
	freeParams=mCommonParameters->namesForPlainValues();
	std::vector<double> initial_params=mCommonParameters->plainValues();
	int numpars=freeParams.size();
	
	//make new columns for these parameters
	for(int i=0; i<numpars; i++){
		mDataTable->addColumn(freeParams[i]);
	}
	
	std::string evalkey="Evaluation";
	mDataTable->addColumn(evalkey);
	double * evalptr=mDataTable->portForColumn(evalkey);
	
	for(int i=0; i<numpars; i++){
		double * port=mDataTable->portForColumn(freeParams[i]);
		if(port){
			paramPorts.push_back(port);
		}
		else{
			printf("[Error]: Could not find data field for parameter %s.\n",freeParams[i].c_str());
		}
	}
	
	//fill parameters with their initial values
	mDataTable->rewind();
	do{
		for(int i=0; i<paramPorts.size(); i++){
			//update parameter i in this row
			*paramPorts[i]=initial_params[i];
		}
	}while(mDataTable->stepRow()!=-1);	
	
	int actstartrow=0;
	int actendrow=-1;
		
	printf("*** Sequential parameter calibration procedure ***\n");
	
	mDataTable->rewind();
	
	//state container
	mModelState.clear();
	
	//iterate over intervals
	while(actstartrow<mDataTable->numRows()){
		
		//search for the next start
		double flagvalue=mDataTable->valueForColumn(eventflagfield, actstartrow);
		actendrow=-1;
		for(int i=actstartrow+1; i<mDataTable->numRows(); i++){
			double next_val=mDataTable->valueForColumn(eventflagfield, i);
			if(flagvalue==0.0 && next_val!=0.0){
				actendrow=i-1;
				break;
			}
			flagvalue=next_val;
		}
		if(actendrow==-1){
			actendrow=mDataTable->numRows()-1;	//end of data table
		}
		if(actendrow<=actstartrow){
			break;
		}
		actstartrow=(actstartrow>0?actstartrow-1:0);
		
		//do that interval
		printf("Optimizing between row #%d and %d...\t",actstartrow,actendrow);
		fflush(stdout);
		
		//set initial values for the parameters
		mCommonParameters->setPlainValues(initial_params);
				
		mEvaluateStartRow=actstartrow;
		mEvaluateEndRow=actendrow;
		calibrate();	
		
		//get best parameter guess
		std::vector<double> best_pars=mCommonParameters->plainValues();
		mEvaluateStartRow=actstartrow;
		mEvaluateEndRow=actendrow;
		double best_guess=evaluate(best_pars);
		printf("Best performace: %lg\n",best_guess);
		
		mModelState=mSolver->modelState();	//at the end of run
		
		//update data table rows with the best guess (previously it was in the evaluate routine)
		if(mDataTable && paramPorts.size()){
			for(int i=mEvaluateStartRow+1; i<=mEvaluateEndRow; i++){
				mDataTable->setRow(i);
				for(int r=0; r<numpars; r++){
					*paramPorts[r]=best_pars[r];	//update by the outlet
				}
				*evalptr=best_guess;
			}
			mDataTable->rewind();	//force commit
		}
		
		//save a backup of the results
		mDataTable->writeToFile("._temporary_results.txt");
		
		//step on the next
		actstartrow=actendrow+1;
	}
	mEvaluateStartRow=-1;
	mEvaluateEndRow=-1;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::sequentialCalibrateInputs(std::string eventflagfield, std::string inputfield, iWQRandomGenerator * inputprior)
{
	//check if eventflagfield points to a valid column
	if(!mDataTable->portForColumn(eventflagfield)){
		printf("[Error]: There is no data column with name: %s\n",eventflagfield.c_str());
	}
	
	if(!mDataTable->portForColumn(inputfield)){
		printf("[Error]: There is no data column with name: %s\n",inputfield.c_str());
	}
	
	//need to have columns for each adjustable parameter in the data table
	std::vector<double *> paramPorts;
	
	std::string evalkey="Evaluation";
	mDataTable->addColumn(evalkey);
	double * evalptr=mDataTable->portForColumn(evalkey);
	
	mRainColPtr=mDataTable->portForColumn(inputfield);
	
	int actstartrow=0;
	int actendrow=-1;
		
	printf("*** Sequential input adjustment procedure ***\n");
	
	mDataTable->rewind();
	
	//state container
	mModelState.clear();
	
	//iterate over intervals
	while(actstartrow<mDataTable->numRows()){
		
		//search for the next start
		double flagvalue=mDataTable->valueForColumn(eventflagfield, actstartrow);
		actendrow=-1;
		for(int i=actstartrow+1; i<mDataTable->numRows(); i++){
			double next_val=mDataTable->valueForColumn(eventflagfield, i);
			if(flagvalue==0.0 && next_val!=0.0){
				actendrow=i-1;
				break;
			}
			flagvalue=next_val;
		}
		if(actendrow==-1){
			actendrow=mDataTable->numRows()-1;	//end of data table
		}
		if(actendrow<=actstartrow){
			break;
		}
		actstartrow=(actstartrow>0?actstartrow-1:0);
		
		//do that interval
		printf("Optimizing between row #%d and %d...\t",actstartrow,actendrow);
		fflush(stdout);
		
		//set initial values for the parameters
		mCommonParameters -> clearAllParams();
		int num_rainy_days = actendrow-actstartrow;
		//create dummy parameters for this interval
		for(int x=0; x<num_rainy_days; x++){
			std::stringstream s;
			s << "R" << x;
			std::string parname = s.str();
			mCommonParameters->initParam(parname, 1.0);
			iWQLimits lim;
			lim.min=0.0;
			lim.max=200.0;
			mCommonParameters->setLimitsForParam(lim, parname, "");
			if(inputprior){
				mCommonParameters->linkDistributionToParam(inputprior, parname);
			}
		}
				
		mEvaluateStartRow=actstartrow;
		mEvaluateEndRow=actendrow;
		
		calibrate();	
		
		//get best parameter guess
		std::vector<double> best_pars=mCommonParameters->plainValues();
		mEvaluateStartRow=actstartrow;
		mEvaluateEndRow=actendrow;
		double best_guess=evaluate(best_pars);
		printf("Best performace: %lg\n",best_guess);
		
		mModelState=mSolver->modelState();	//at the end of run
		
		//save a backup of the results
		mDataTable->writeToFile("._temporary_results.txt");
		
		//step on the next
		actstartrow=actendrow+1;
	}
	mEvaluateStartRow=-1;
	mEvaluateEndRow=-1;
	mRainColPtr=NULL;
}

//-----------------------------------------------------------------------------------

void iWQEvaluator::calibrate()
{
	//Heavy reuse of asa047 code in this function
	int i, j;
	int icount;
	int ifault;
	int kcount;
	int kround;
	int konvge;
	int n;
	int numres;
	double reqmin;
	double * start;
	double * step;
	double * xmin;
	double ynewlo;
	
	n = mCommonParameters->numberOfParams();
	
	start = new double[n];
	step = new double[n];
	xmin = new double[n];
	
	reqmin = NMSTolerance;
	
	konvge = 10;
	kcount = NMSMaxNumRounds;	//maximal number of rounds
	kround = 100;	//number of iterations in each round
	double stepfactor=1.0;
	double lasty;
	
	std::vector<double> parvals;
	std::vector<std::string> parnames=mCommonParameters->namesForPlainValues();
	
	printWarnings=false;
	
	//preparatory phase with PSO
	if(PSOActive){
		printf("Particle Swarm Optimization...\n");
		parvals=mCommonParameters->plainValues();
		
		iWQBoundsList bounds;
		for(i=0; i<parvals.size(); i++){
			if(mCommonParameters->hasLimitsForParam(parnames[i])){
				//we have bounds
				iWQLimits lim=mCommonParameters->limitsForParam(parnames[i]);
				bounds.add(lim.min,lim.max);
			}
			else{
				//no defined bounds
				bounds.add(0,(parvals[i]!=0.0?10*parvals[i]:0.0));		//search from 0 to 10*parvals[i]
			}
		}
		parvals=iWQParticleSwarmOptimize(this,bounds,PSOSwarmSize,PSOMaxNumRounds,PSOMaxIdleRounds);
		mCommonParameters->setPlainValues(parvals);
		printf("Ready\n");
	}
	
	if(NMSActive){
		printf("Nelder-Mead Simplex Optimization...\n");
		for(j=0; j<kcount; j++){
			parvals=mCommonParameters->plainValues();
			for(i=0; i<n; i++){
				start[i]=parvals[i];
			}
			for(i=0; i<n; i++){
				step[i]=stepfactor * (parvals[i]!=0.0?parvals[i]:0.1) / 5.0;	//20% variation + prevent stucking if a parameter is 0
			}
		
			NelderMead( n, start, xmin, &ynewlo, reqmin, step, konvge, kround, &icount, &numres, &ifault );
		
			for(i=0; i<n; i++){
				parvals[i]=xmin[i];
			}
			stepfactor*=0.99;
			mCommonParameters->setPlainValues(parvals);
			printf("[%0.6g]\n",ynewlo);
		
			//save parameters to a temporary file
			FILE * tempfile=fopen("_calibration_progress.tmp","a");
			if(tempfile){
				time_t t=time(0);
				fprintf(tempfile,"#BEGIN RECORD\n#time=%s\n#creator=simplex\n#iteration=%d\n",ctime(&t),j);
				for(i=0; i<parvals.size(); i++){
					fprintf(tempfile,"\t%s: %g\n",parnames[i].c_str(),parvals[i]);
				}
				fprintf(tempfile,"#eval=[%g]\n#END RECORD\n",ynewlo);
				fclose(tempfile);
			}
			//mCommonParameters->saveToFile("_calibration_progress.tmp");
		
			if(j>0){
				if(fabs(ynewlo-lasty)<reqmin || ynewlo>lasty){
					//printf("Iteration stuck.\n");
					break;
				}
			}
			lasty=ynewlo;
		}
	}
	printWarnings=true;
	delete [] start;
	delete [] step;
	delete [] xmin;
}

//#######################################################################################

void iWQEvaluator::NelderMead(int n, 
							  double start[], 
							  double xmin[], 
							  double *ynewlo, 
							  double reqmin, 
							  double step[], 
							  int konvge, 
							  int kcount, 
							  int *icount, 
							  int *numres, 
							  int *ifault)
{	
	//  Nelder-Mead algorithm from asa047.C     
	//  Author: C++ version by John Burkardt	
	double ccoeff = 0.5;
	double del;
	double dn;
	double dnn;
	double ecoeff = 2.0;
	double eps = 0.001;
	int i;
	int ihi;
	int ilo;
	int j;
	int jcount;
	int l;
	int nn;
	double *p;
	double *p2star;
	double *pbar;
	double *pstar;
	double rcoeff = 1.0;
	double rq;
	double x;
	double *y;
	double y2star;
	double ylo;
	double ystar;
	double z;
	//
	//  Check the input parameters.
	//
	if ( reqmin <= 0.0 )
	{
		*ifault = 1;
		return;
	}
	
	if ( n < 1 )
	{
		*ifault = 1;
		return;
	}
	
	if ( konvge < 1 )
	{
		*ifault = 1;
		return;
	}
	
	p = new double[n*(n+1)];
	pstar = new double[n];
	p2star = new double[n];
	pbar = new double[n];
	y = new double[n+1];
	
	*icount = 0;
	*numres = 0;
	
	jcount = konvge; 
	dn = ( double ) ( n );
	nn = n + 1;
	dnn = ( double ) ( nn );
	del = 1.0;
	rq = reqmin * dn;
	//
	//  Initial or restarted loop.
	//
	for ( ; ; )
	{
		for ( i = 0; i < n; i++ )
		{ 
			p[i+n*n] = start[i];
		}
		y[n] = evaluate( start, n );
		*icount = *icount + 1;
		
		for ( j = 0; j < n; j++ )
		{
			x = start[j];
			start[j] = start[j] + step[j] * del;
			for ( i = 0; i < n; i++ )
			{
				p[i+j*n] = start[i];
			}
			y[j] = evaluate( start, n );
			*icount = *icount + 1;
			start[j] = x;
		}
		//                    
		//  The simplex construction is complete.
		//                    
		//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
		//  the vertex of the simplex to be replaced.
		//                
		ylo = y[0];
		ilo = 0;
		
		for ( i = 1; i < nn; i++ )
		{
			if ( y[i] < ylo )
			{
				ylo = y[i];
				ilo = i;
			}
		}
		//
		//  Inner loop.
		//
		for ( ; ; )
		{
			if ( kcount <= *icount )
			{
				break;
			}
			*ynewlo = y[0];
			ihi = 0;
			
			for ( i = 1; i < nn; i++ )
			{
				if ( *ynewlo < y[i] )
				{
					*ynewlo = y[i];
					ihi = i;
				}
			}
			//
			//  Calculate PBAR, the centroid of the simplex vertices
			//  excepting the vertex with Y value YNEWLO.
			//
			for ( i = 0; i < n; i++ )
			{
				z = 0.0;
				for ( j = 0; j < nn; j++ )
				{ 
					z = z + p[i+j*n];
				}
				z = z - p[i+ihi*n];  
				pbar[i] = z / dn;
			}
			//
			//  Reflection through the centroid.
			//
			for ( i = 0; i < n; i++ )
			{
				pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
			}
			ystar = evaluate( pstar, n );
			*icount = *icount + 1;
			//
			//  Successful reflection, so extension.
			//
			if ( ystar < ylo )
			{
				for ( i = 0; i < n; i++ )
				{
					p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
				}
				y2star = evaluate( p2star, n );
				*icount = *icount + 1;
				//
				//  Check extension.
				//
				if ( ystar < y2star )
				{
					for ( i = 0; i < n; i++ )
					{
						p[i+ihi*n] = pstar[i];
					}
					y[ihi] = ystar;
				}
				//
				//  Retain extension or contraction.
				//
				else
				{
					for ( i = 0; i < n; i++ )
					{
						p[i+ihi*n] = p2star[i];
					}
					y[ihi] = y2star;
				}
			}
			//
			//  No extension.
			//
			else
			{
				l = 0;
				for ( i = 0; i < nn; i++ )
				{
					if ( ystar < y[i] )
					{
						l = l + 1;
					}
				}
				
				if ( 1 < l )
				{
					for ( i = 0; i < n; i++ )
					{
						p[i+ihi*n] = pstar[i];
					}
					y[ihi] = ystar;
				}
				//
				//  Contraction on the Y(IHI) side of the centroid.
				//
				else if ( l == 0 )
				{
					for ( i = 0; i < n; i++ )
					{
						p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
					}
					y2star = evaluate( p2star, n );
					*icount = *icount + 1;
					//
					//  Contract the whole simplex.
					//
					if ( y[ihi] < y2star )
					{
						for ( j = 0; j < nn; j++ )
						{
							for ( i = 0; i < n; i++ )
							{
								p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
								xmin[i] = p[i+j*n];
							}
							y[j] = evaluate( xmin, n );
							*icount = *icount + 1;
						}
						ylo = y[0];
						ilo = 0;
						
						for ( i = 1; i < nn; i++ )
						{
							if ( y[i] < ylo )
							{
								ylo = y[i];
								ilo = i;
							}
						}
						continue;
					}
					//
					//  Retain contraction.
					//
					else
					{
						for ( i = 0; i < n; i++ )
						{
							p[i+ihi*n] = p2star[i];
						}
						y[ihi] = y2star;
					}
				}
				//
				//  Contraction on the reflection side of the centroid.
				//
				else if ( l == 1 )
				{
					for ( i = 0; i < n; i++ )
					{
						p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
					}
					y2star = evaluate( p2star, n );
					*icount = *icount + 1;
					//
					//  Retain reflection?
					//
					if ( y2star <= ystar )
					{
						for ( i = 0; i < n; i++ )
						{
							p[i+ihi*n] = p2star[i];
						}
						y[ihi] = y2star;
					}
					else
					{
						for ( i = 0; i < n; i++ )
						{
							p[i+ihi*n] = pstar[i];
						}
						y[ihi] = ystar;
					}
				}
			}
			//
			//  Check if YLO improved.
			//
			if ( y[ihi] < ylo )
			{
				ylo = y[ihi];
				ilo = ihi;
			}
			jcount = jcount - 1;
			
			if ( 0 < jcount )
			{
				continue;
			}
			//
			//  Check to see if minimum reached.
			//
			if ( *icount <= kcount )
			{
				jcount = konvge;
				
				z = 0.0;
				for ( i = 0; i < nn; i++ )
				{
					z = z + y[i];
				}
				x = z / dnn;
				
				z = 0.0;
				for ( i = 0; i < nn; i++ )
				{
					z = z + pow( y[i] - x, 2 );
				}
				
				if ( z <= rq )
				{
					break;
				}
			}
		}
		//
		//  Factorial tests to check that YNEWLO is a local minimum.
		//
		for ( i = 0; i < n; i++ )
		{
			xmin[i] = p[i+ilo*n];
		}
		*ynewlo = y[ilo];
		
		if ( kcount < *icount )
		{
			*ifault = 2;
			break;
		}
		
		*ifault = 0;
		
		for ( i = 0; i < n; i++ )
		{
			del = step[i] * eps;
			xmin[i] = xmin[i] + del;
			z = evaluate( xmin, n );
			*icount = *icount + 1;
			if ( z < *ynewlo )
			{
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] - del - del;
			z = evaluate( xmin, n );
			*icount = *icount + 1;
			if ( z < *ynewlo )
			{
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] + del;
		}
		
		if ( *ifault == 0 )
		{
			break;
		}
		//
		//  Restart the procedure.
		//
		for ( i = 0; i < n; i++ )
		{
			start[i] = xmin[i];
		}
		del = eps;
		*numres = *numres + 1;
	}
	delete [] p;
	delete [] pstar;
	delete [] p2star;
	delete [] pbar;
	delete [] y;
	
	return;
}
