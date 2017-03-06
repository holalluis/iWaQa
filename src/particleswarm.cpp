/*
 *  particleswarm.cpp
 *  Particle-swarm optimiser
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/OPTIMISE
 *
 */ 
 
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "particleswarm.h"
#include "evaluator.h"
#include "model.h"

//----------------------------------------------------------------------------

void iWQBoundsList::add(double min, double max)
{ 
	push_back(iWQMakeBounds((min<max?min:max), (max>=min?max:min))); 
}

//----------------------------------------------------------------------------

iWQBounds iWQMakeBounds(double _min, double _max)
{
	iWQBounds b;
	b.min=_min;
	b.max=_max;
	return b;
}

//----------------------------------------------------------------------------

//auxiliary functions
template<class T>
T ** make2Darray(int dim1, int dim2)
{
	T ** result;
	result=new T * [dim1];
	for(int i=0; i<dim1; i++){
		result[i]=new T [dim2];
		memset(result[i],0,dim2*sizeof(T));
	}
	return result;
}

//----------------------------------------------------------------------------

template<class T>
void delete2Darray(T ** array, int dim1, int dim2)
{
	if(!array){
		return; 
	}
	for(int i=0; i<dim1; i++){
		delete [] array[i];
	}
	delete [] array;
}

//----------------------------------------------------------------------------

template<class T>
T * make1Darray(int dim1)
{
	T * result;
	result=new T [dim1];
	memset(result,0,dim1*sizeof(T));
	return result;
}

//----------------------------------------------------------------------------

double Rnd()
{
	return rand()/((double)RAND_MAX + 1.0);
}

//----------------------------------------------------------------------------

//optimizer implementation (translated from VISUAL BASIC implementation at http://read.pudn.com/downloads137/sourcecode/math/587436/FRMSWARM.FRM__.htm)

std::vector<double> iWQParticleSwarmOptimize(iWQEvaluator * evaluator, iWQBoundsList bounds, int populationsize, int maxiterations, int idlerunlength)
{
	//declarations
	int iPOPSIZE; 		// Population size
	int iDIMENSIONS; 	// Number of dimensions
	double fINITWT;		// Initial inertia weight
	double fMAXVEL; 	// Maximum velocity allowed
	int nMAXITER;		// Maximum number of iterations
	double ** fPos;		// Position for each particle
	double ** fTempPos; // Temporary updated position for each particle
	double ** fVel;		// Velocity for each particle
	double * fDumVel;	// Dummy velocity vector 1-dim array
	double ** fBestPos;	// Best previous position for each particle
	iWQBoundsList fBounds;	//bounds for each dimension
	double fInerWt; 	// Inertia weight
	double * fErrVal; 	// Function error value calculated by eval()
	double * fPbestVal;	// Best error value over time for each particle
	double fERRCUTOFF; 	// Error value at which system stops
	int iUSEBETTER;		// Usebetter is kind of a momentum
	int * iBetter; 		// This gets set in program
	int iLOCAL;			// Neighborhood size specified in run file
	int iHOODSIZE;		// Neighborhood size used in program
	int * iNeighbor;	// Popindexes of neighbors resulting from iLOCAL
	
	printf("Running particle swarm optimization test.\n");
	
	if(!evaluator){
		printf("[Error]: no evaluator specified for PSO optimization.\n");
		return std::vector<double> ();
	}
	
	//init things
	iDIMENSIONS=bounds.size();
	
	//get parameter names
	iWQParameterManager * parmanager = evaluator->parameters();
	
	if(!parmanager){
		printf("[Error]: Evaluator does not have parameters for PSO optimization.\n");
		return std::vector<double> ();
	}
	
	std::vector<std::string> parnames = parmanager->namesForPlainValues();
	
	//get initial parameter values from parmanager
	std::vector<double> parvals = parmanager->plainValues();
		
	if(iDIMENSIONS>parnames.size()){
		printf("[Error]: The count of parameter names does not match PSO bounds dimension.\n");
		return std::vector<double> ();
	}
	
	//other parameters
	iPOPSIZE = populationsize; 	
	fERRCUTOFF = idlerunlength;			//1E-6
	fMAXVEL = 0.1;				//needs scaling in each dimension	//10
	//fMaxPos = 500;				//100
	nMAXITER = maxiterations;			//2000
	iUSEBETTER = 0;				//0
	fINITWT = 0.9;				//0.9
	iLOCAL = 2 * iDIMENSIONS;				//2
	
	iHOODSIZE = (iLOCAL>0)?((iLOCAL / 2) * 2):iPOPSIZE;
	if(iHOODSIZE>iPOPSIZE){
		iHOODSIZE=iPOPSIZE;
	}
	
	//furnish storages
	fPos = make2Darray<double>(iPOPSIZE, iDIMENSIONS);
	fTempPos = make2Darray<double>(iPOPSIZE, iDIMENSIONS);
	fVel = make2Darray<double>(iPOPSIZE, iDIMENSIONS);
	fBestPos = make2Darray<double>(iPOPSIZE, iDIMENSIONS);
	fDumVel = make1Darray<double>(iDIMENSIONS);
	fErrVal = make1Darray<double>(iPOPSIZE);
	fPbestVal = make1Darray<double>(iPOPSIZE);
	iBetter = make1Darray<int>(iPOPSIZE);
	iNeighbor = make1Darray<int>(iHOODSIZE+1); //! -iPOPSIZE to iPOPSIZE
  
	//init randomization
	srand(time(0));
	
	int iPopindex; 		// Index for population
	int iDimindex;		// Index for dimensions
	int nIter;			// Number of iterations
	int iHOODINDEX;		// Neighborhood index (offset from particle)
	int iGbest; 		// Index for global best particle
	int iLbest; 		// Index for local (neighborhood) best particle
	double previousbest=0.0;
	int samebestcount=0;
	double modelpos [iDIMENSIONS];
	
	fBounds=bounds;
		
	// Randomize the positions and velocities for entire population
	for(iPopindex=0; iPopindex<iPOPSIZE; iPopindex++){
		for(iDimindex = 0;  iDimindex<iDIMENSIONS; iDimindex++){
			if(iPopindex==0){
				fPos[iPopindex][iDimindex] = (parvals[iDimindex] - fBounds[iDimindex].min) / (fBounds[iDimindex].max-fBounds[iDimindex].min);	//the 1st particle is positioned in the original parameter value
			}
			else{
				fPos[iPopindex][iDimindex] = Rnd(); // * fMaxPos;
            }
			fBestPos[iPopindex][iDimindex] = fPos[iPopindex][iDimindex];
            fVel[iPopindex][iDimindex] = Rnd() * fMAXVEL;
            if(Rnd()>0.5){
				fVel[iPopindex][iDimindex] *= -1.0; 
			}
		}
	}

	// Main swarm loop here
    for(nIter=0; nIter<nMAXITER; nIter++){
    
        // Update inertia weight; linear from fINITWT to 0.4
        fInerWt = ((fINITWT - 0.4) * (nMAXITER - nIter) / (double)nMAXITER) + 0.4;
        
		for(iPopindex = 0;  iPopindex<iPOPSIZE; iPopindex++){     			// MAIN main loop starts here
            // Setup dummy velocity vector for current population member
            for(iDimindex = 0; iDimindex<iDIMENSIONS; iDimindex++){
				fDumVel[iDimindex] = fVel[iPopindex][iDimindex];
            }
            
            iBetter[iPopindex] = 0;             							// Set to 0 unless new Pbest achieved
            
			//translate to model space 
			for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){
				modelpos[iDimindex]=fBounds[iDimindex].min+fPos[iPopindex][iDimindex]*(fBounds[iDimindex].max-fBounds[iDimindex].min);
			}
			
            fErrVal[iPopindex]=evaluator->evaluate(modelpos, iDIMENSIONS);    			// evaluates f function: fErrVal(iPopindex)
            
            if(nIter==0){
                fPbestVal[iPopindex] = fErrVal[iPopindex];
                iGbest = 0;
            }
            
            if(fErrVal[iPopindex] < fPbestVal[iPopindex]){   				//If new Pbest
                fPbestVal[iPopindex] = fErrVal[iPopindex];
                
                for(iDimindex = 0; iDimindex<iDIMENSIONS; iDimindex++){    	// Reset Pbest location vector
                    fBestPos[iPopindex][iDimindex] = fPos[iPopindex][iDimindex];
                }
                
                if(fPbestVal[iPopindex] < fPbestVal[iGbest]){
                    iGbest = iPopindex;
                }
                
                if(iUSEBETTER==1){
					iBetter[iPopindex] = 1;
				}
            }                                     							// End new Pbest condition
		}								          							// end MAIN main loop for gold gbest only
		
		for(iPopindex = 0; iPopindex<iPOPSIZE; iPopindex++){         		// update velocity & position
      
            // Does neighborhood calculation of iLbest
            if(iLOCAL > 0){
				//TODO: check this
                for(iHOODINDEX = 0; iHOODINDEX<=iHOODSIZE; iHOODINDEX++){
                    iNeighbor[iHOODINDEX] = iPopindex - (iHOODSIZE / 2) + iHOODINDEX;
                    // for iPopindex = 1,goes from 0 to 2 for iHOODSIZE of 2
                    //                       from -1 to 3 for iHOODSIZE of 4
                    
					// Now wrap the ends of the array
                    if(iNeighbor[iHOODINDEX] < 0){
						iNeighbor[iHOODINDEX] += iPOPSIZE;
					}
                    if(iNeighbor[iHOODINDEX] >= iPOPSIZE){
						iNeighbor[iHOODINDEX] -= iPOPSIZE;
					}
                    // Start with iNeighbor[0] as iLbest and try to beat it
                    if(iHOODINDEX == 0){
						iLbest = iNeighbor[0];
					}
                    if(fPbestVal[iNeighbor[iHOODINDEX]] < fPbestVal[iLbest]){
						iLbest = iNeighbor[iHOODINDEX];
					}
                }
            }
            
			if(iLOCAL == 0){
				iLbest = iGbest;
			}
                        
            // Update velocity vector for one particle Russ Reduced version
            for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){        	//fInerWt below
                fVel[iPopindex][iDimindex] = (0.5 + (Rnd() / 2.0)) * fVel[iPopindex][iDimindex] + 2.0 * Rnd() * 
											 (fBestPos[iPopindex][iDimindex] - fPos[iPopindex][iDimindex]) + 
											 2.0 * Rnd() * (fBestPos[iLbest][iDimindex] - fPos[iPopindex][iDimindex]);
                 
                if(fVel[iPopindex][iDimindex] > fMAXVEL){
                    fVel[iPopindex][iDimindex] = fMAXVEL;
				}
                else if(fVel[iPopindex][iDimindex] < -fMAXVEL){
                    fVel[iPopindex][iDimindex] = -fMAXVEL;
                }
            }
            
            //If it's going the right way, keep going
            if(iBetter[iPopindex]== 1){
                for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){
                    fVel[iPopindex][iDimindex] = fDumVel[iDimindex];
                }
            }

            for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){  		// Define new positions for all dimensions
                fPos[iPopindex][iDimindex] = fPos[iPopindex][iDimindex] + fVel[iPopindex][iDimindex];
            }
		}                	//end velocity & position loop
		
		//write out iteration details
		//translate to model space
		for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){
			modelpos[iDimindex]=fBounds[iDimindex].min+fBestPos[iGbest][iDimindex]*(fBounds[iDimindex].max-fBounds[iDimindex].min);
		}
		//write best position
		printf("PSO #%d",nIter);
		//write best value
		printf("\t[%lf]\n",fPbestVal[iGbest]);
		
        // Terminate on sufficiently low error
        if(nIter>0 && previousbest==fPbestVal[iGbest]){
			samebestcount++;
		}
		else{
			samebestcount=0;
		}
		previousbest=fPbestVal[iGbest];
		if((samebestcount >= fERRCUTOFF) || (nIter >= nMAXITER)){
			//translate to model space
			break;				//END PROGRAM
        }
		
		//save parameters to a temporary file
		FILE * tempfile=fopen("_calibration_progress.tmp","a");
		if(tempfile){
			time_t t=time(0);
			fprintf(tempfile,"#BEGIN RECORD\n#time=%s\n#creator=particleswarm\n#iteration=%d\n",ctime(&t),nIter);
			for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){
				fprintf(tempfile,"\t%s: %g\n",parnames[iDimindex].c_str(), modelpos[iDimindex]);
			}
			fprintf(tempfile,"#eval=[%g]\n#END RECORD\n",fPbestVal[iGbest]);
			fclose(tempfile);
		}
	}                       // next nIter, end main loop
	//fclose(ofile);
	std::vector<double> result (iDIMENSIONS);
	for(iDimindex=0; iDimindex<iDIMENSIONS; iDimindex++){
		modelpos[iDimindex]=fBounds[iDimindex].min+fBestPos[iGbest][iDimindex]*(fBounds[iDimindex].max-fBounds[iDimindex].min);
		result[iDimindex]=modelpos[iDimindex];
	}
	
	//delete storages
	delete2Darray<double>(fPos, iPOPSIZE, iDIMENSIONS);
	delete2Darray<double>(fTempPos, iPOPSIZE, iDIMENSIONS);
	delete2Darray<double>(fVel, iPOPSIZE, iDIMENSIONS);
	delete2Darray<double>(fBestPos, iPOPSIZE, iDIMENSIONS);
	delete [] fDumVel;
	delete [] fErrVal;
	delete [] fPbestVal;
	delete [] iBetter;
	delete [] iNeighbor;
	
	return result;
}
