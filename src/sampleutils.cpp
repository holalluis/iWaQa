/*
 *  sampleutils.cpp
 *  Histogram creator from MCMC samples
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/UTILS
 *
 */ 

#include <string>
#include <vector>
#include <stdio.h>

#include "model.h"
#include "datatable.h"
#include "mathutils.cpp"
  
void iWQEstimateDistributions(std::string samplefilename, iWQStrings paramnames, std::string outputfilename, int burn_in_length)
{
	iWQDataTable * datatable=new iWQDataTable (samplefilename);
	
	//check if we have anything in the sample
	if(datatable->numRows()<=0){
		printf("[Error]: Failed to load sample data from file \"%s\".\n",samplefilename.c_str());
		delete datatable;
		return;
	}
	
	int nrows=datatable->numRows();
	int nparams=paramnames.size();
	
	//make flag storage (indicates burn-in and repetitions with "true")
	std::vector<bool> flags;
	flags.assign(nrows, false);
	
	//now check parameter repetitions
	std::vector< std::vector<int> > repflags;
	std::vector< const std::vector<double> * > datacols;
	int maxreplimit=5;
	for(int i=0; i<nparams; i++){
		std::vector<int> v;
		v.assign(nrows, 0);
		repflags.push_back(v);
		std::string act_param=paramnames[i];
		const std::vector<double> * rawdata=datatable->vectorForColumn(act_param);
		datacols.push_back(rawdata);
	}
	for(int p=0; p<nparams; p++){
		double prev_val=datacols[p]->at(0);
		double prev_repflag=0;
		for(int r=1; r<nrows; r++){
			double act_val=datacols[p]->at(r);
			double act_repflag=0;
			if(act_val==prev_val){
				act_repflag=prev_repflag+1;
			}
			repflags[p][r]=(int)act_repflag;
			prev_repflag=act_repflag;
			prev_val=act_val;
		}
	}
	for(int r=0; r<nrows; r++){
		bool rep=true;
		for(int p=0; p<nparams; p++){
			if(repflags[p][r]<maxreplimit){
				rep=false;
			}
		}
		if(rep){
			flags[r]=true;
		}
	}
	
	FILE * ofile=fopen(outputfilename.c_str(), "w");
	if(!ofile){
		printf("[Error]: Failed to open output file \"%s\".\n",outputfilename.c_str());
		delete datatable;
		return;
	}
	
	//now iterate through all parameters
	for(int i=0; i<nparams; i++){
		
		std::string act_param=paramnames[i];
		
		std::vector<double> * data=new std::vector<double>;
		for(int r=burn_in_length; r<nrows; r++){
			if(!flags[r]){
				data->push_back(datacols[i]->at(r));
			}
		}
		
		int ndata=data->size();
		
		if(data && ndata>0){
			//get min and max
			double Min=min(data);
			double Max=max(data);
						
			fprintf(ofile,"\nHISTOGRAM OF %s\n",act_param.c_str());
			fprintf(ofile,"Bin_start\tBin_end\tBin_middle\tCount\tProportion\n");
				
			if(Min<DBL_MAX && Max>-DBL_MAX){
				//confined interval
				int numinvaliddata=0;	//extra bin for INFs and NaNs
				std::vector<int> counts;
				
				if(Max>Min){
					int steps=(ndata>500)?50:ndata/10;
					double stepsize=(Max-Min)/(double)steps;
					counts.assign(steps+1,0);
					for(int i=0; i<ndata; i++){
						double val=data->at(i);
						if(!isnan(val) && !isinf(val)){
							int binindex=(int)((val-Min)/stepsize);
							counts[binindex]++;
						}
						else{
							numinvaliddata++;
						}						
					}
					
					for(int i=0; i<counts.size(); i++){
						fprintf(ofile,"%g\t%g\t%g\t%d\t%g\n",Min+(double)i*stepsize,Min+(i+1.0)*stepsize,Min+(i+0.5)*stepsize, counts[i], counts[i]/(double)ndata);
					}
				}
				else{
					//a single value
					counts.assign(1,0);
					for(int i=0; i<ndata; i++){
						double val=data->at(i);
						if(!isnan(val) && !isinf(val)){
							counts[0]++;
						}
						else{
							numinvaliddata++;
						}						
					}
					fprintf(ofile,"%g\t%g\t%g\t%d\t%g\n",Min,Min,Min, counts[0], counts[0]/(double)ndata);
				}
				//apend the invalid data row
				fprintf(ofile,"#INVALID_DATA\t\t\t%d\t%g\n",numinvaliddata, numinvaliddata/(double)ndata);
			}
			else{
				fprintf(ofile,"#INVALID_DATA\t\t\t%d\t%g\n",ndata, 1.0);
			}
			
			delete data;
		}		
	}
	fclose(ofile);
	delete datatable;
}

