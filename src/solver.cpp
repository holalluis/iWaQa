/*
 *  solver.cpp
 *  General model solver
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MODEL
 *
 */ 

#include <stdio.h>
 
#include "solver.h"
#include "datatable.h"

//--------------------------------------------------------------------------------------------------

#pragma mark Link

void iWQLink::init()
{
	srcmod=NULL;
	destmod=NULL; 
	srcptr=NULL;
	destptr=NULL;
		
	fixed_proportion=false;
	proportion=1.0;
	keyed_proportion=false;
	prop_numerator=NULL;
	prop_denominator=NULL;
}

//--------------------------------------------------------------------------------------------------

void iWQLink::establish(iWQModel * src, std::string srcname, iWQModel * dest, std::string destname)
{ 
	init();
	if(!src || !dest){
		return;
	}
	srcmod=src; 
	srcptr=src->routlet(srcname);
	destmod=dest; 
	destptr=destmod->rwoutlet(destname);
}

//--------------------------------------------------------------------------------------------------

void iWQLink::establish(const double * src, iWQModel * dest, std::string destname)
{
	init();
	if(!src || !dest){
		return;
	}
	srcptr=src;
	destmod=dest;
	destptr=destmod->rwoutlet(destname);
}

//--------------------------------------------------------------------------------------------------

void iWQLink::establish(iWQModel * src, std::string srcname, double * dest)
{
	init();
	if(!src || !dest){
		return;
	}
	srcmod=src; 
	srcptr=src->routlet(srcname);
	destptr=dest;
}

//--------------------------------------------------------------------------------------------------

void iWQLink::establish(const double * src, double * dest)
{
	init();
	if(!src || !dest){
		return;
	}
	srcptr=src;
	destptr=dest;
}

//--------------------------------------------------------------------------------------------------

void iWQLink::zerodest()
{
	if(destptr){
		*destptr=0.0;
	}
	
	//refresh dynamic proportion
	if(keyed_proportion){
		if(*prop_denominator!=0.0){
			proportion=(*prop_numerator)/(*prop_denominator);
		}
		else{
			proportion=1.0;
		}
	}
}

//--------------------------------------------------------------------------------------------------

void iWQLink::linkadd()
{
	if(destptr && srcptr){ 
		*destptr += proportion * (*srcptr); 	
	}
}

//--------------------------------------------------------------------------------------------------

iWQModel * iWQLink::dependsOn()
{ 
	return srcmod; 
}

//--------------------------------------------------------------------------------------------------

iWQModel * iWQLink::subject()
{ 
	return destmod; 
}

//--------------------------------------------------------------------------------------------------

void iWQLink::print()
{
	printf("[iWQLink from 0x%p (0x%p) to 0x%p (0x%p)]",srcmod,srcptr,destmod,destptr);
}

//--------------------------------------------------------------------------------------------------

void iWQLink::setFixedProportion(double prop)
{
	fixed_proportion=true;
	keyed_proportion=false;
	proportion=prop;
	prop_numerator=NULL;
	prop_denominator=NULL;
}

//--------------------------------------------------------------------------------------------------
		
void iWQLink::setKeyedProportion(std::string key)
{
	fixed_proportion=false;
	keyed_proportion=false;
	proportion=1.0;
	prop_numerator=NULL;
	prop_denominator=NULL;
	//try to get pointers from src and dest
	if(srcmod && destmod){
		keyed_proportion=true;
		prop_numerator=destmod->routlet(key);
		prop_denominator=srcmod->routlet(key);
	}
	if(prop_numerator && prop_denominator){
		keyed_proportion=true;
	}
	else{
		printf("[Error]: Could not make link proportional to %s.\n",key.c_str());
	}
}	

//--------------------------------------------------------------------------------------------------

bool iWQLink::operator==(const iWQLink & alink) const
{
	return (
		srcmod==alink.srcmod &&
		destmod==alink.destmod && 
		srcptr==alink.srcptr && 
		destptr==alink.destptr && 	
		fixed_proportion==alink.fixed_proportion && 
		proportion==alink.proportion &&
		keyed_proportion==alink.keyed_proportion && 
		prop_numerator==alink.prop_numerator &&
		prop_denominator==alink.prop_denominator
	);
}

//--------------------------------------------------------------------------------------------------

bool iWQLink::operator!=(const iWQLink & alink) const
{
	return !(*this == alink);
}
		
//##################################################################################################

#pragma mark Solver

iWQSolver::iWQSolver(iWQLinkSet links, iWQLinkSet outputlinks)
{
	mTreeError=false;
	
	if(links.size()==0 && outputlinks.size()==0){
		//nothing to do
		return;
	}
	
	//default precision 
	mHmin=1.0/1440.0;	//minute resolution on a daily scale
	mEps=0.001;
	
	std::vector<iWQModel *> modelbuf;
		
	//get model list in mModels
	for(int i=0; i<links.size(); i++){
		//look for units in destinations
		iWQModel * act_model=links[i].subject();
		if(act_model){
			bool alreadyinside=false;
			for(int j=0; j<modelbuf.size(); j++){
				if(modelbuf[j]==act_model){
					alreadyinside=true;
					break;
				}
			}
			if(!alreadyinside){
				modelbuf.push_back(act_model);
			}
		}
		//look for those, which are only sources
		act_model=links[i].dependsOn();
		if(act_model){
			bool alreadyinside=false;
			for(int j=0; j<modelbuf.size(); j++){
				if(modelbuf[j]==act_model){
					alreadyinside=true;
					break;
				}
			}
			if(!alreadyinside){
				modelbuf.push_back(act_model);
			}
		}
	}
	
	//check if a unit is only related to output (has no input or other connections)
	for(int i=0; i<outputlinks.size(); i++){
		//look for units in destinations
		iWQModel * act_model=outputlinks[i].dependsOn();
		if(act_model){
			bool alreadyinside=false;
			for(int j=0; j<modelbuf.size(); j++){
				if(modelbuf[j]==act_model){
					alreadyinside=true;
					break;
				}
			}
			if(!alreadyinside){
				modelbuf.push_back(act_model);
			}
		}
	}
	
	//DEBUG
	//printf("Found %zd model units:\n",modelbuf.size());
	for(int i=0; i<modelbuf.size(); i++){
		iWQModel * mmodel=modelbuf[i];
		//printf("\t#%d\t%s (%s)\n",i+1,mmodel->modelId().c_str(),mmodel->modelType().c_str());
	}
	//END DEBUG
	
	//simple assignment for input links
	mLinks=links;
	
	//also for output links
	mExportLinks=outputlinks;
	
	//select inter-model links
	mInterLinks.clear();
	for(int i=0; i<mLinks.size(); i++){
		const iWQModel * source=mLinks[i].srcmod;
		const iWQModel * dest=mLinks[i].destmod;
		if(source && dest){
			mInterLinks.push_back(mLinks[i]);
		}	
	} 
	
	//get tree root
	//const iWQModel * tree_root=NULL;
	//int tree_root_index=-1;
	std::vector<int> tree_roots;
	tree_roots.clear();
	
	std::vector<int> numDependentNeighbours;
	numDependentNeighbours.assign(modelbuf.size(),0);
	for(int i=0; i<modelbuf.size(); i++){
		const iWQModel * act_model=modelbuf[i];
		
		for(int j=0; j<mLinks.size(); j++){
			if(mLinks[j].dependsOn()==act_model){
				numDependentNeighbours[i]++;
			}
		}
	}
	
	for(int i=0; i<numDependentNeighbours.size(); i++){
		if(numDependentNeighbours[i]==0){
			//tree_root=modelbuf[i];
			//tree_root_index=i;
			tree_roots.push_back(i);
		}
	}
	
	if(tree_roots.size()==0/*tree_root_index==-1*/){
		printf("[Error]: Could not find the root of the model tree.\n");
		mTreeError=true;
		return;
	}
	
	//DEBUG
	else{
		if(tree_roots.size()==1){
			//printf("The root element is:\n");
		}
		else{
			//printf("The root elements are:\n");
		}
		for(int i=0; i<tree_roots.size(); i++){
			const iWQModel * mmodel=modelbuf[tree_roots[i]];
			//printf("#%d %s (%s)\n",tree_roots[i]+1,mmodel->modelId().c_str(),mmodel->modelType().c_str());
		}
	}
	//END DEBUG
	
	//try to find tree layers for each tree
	std::vector<int> layerIndex;
	layerIndex.assign(modelbuf.size(),-1);
	int max_layer_index=-1;
	
	for(int tr=0; tr<tree_roots.size(); tr++){
		layerIndex[tree_roots[tr]]=0;
		
		//without recursion: multiple parse runs
		int act_index=0;
		int numberofdescendants;
		do{
			numberofdescendants=0;
			for(int i=0; i<modelbuf.size(); i++){
				if(layerIndex[i]==act_index){
					//the layerIndex of this unit matches the currently seeked value
					for(int j=0; j<mLinks.size(); j++){		//examine each link 
						if(mLinks[j].subject()==modelbuf[i]){
							//found a link from the current unit to...
							const iWQModel * descendant=mLinks[j].dependsOn();	//descendant should be at a higher index later
							for(int k=0; k<modelbuf.size(); k++){
								if(modelbuf[k]==descendant){
									if(layerIndex[k]!=-1){
										//we have been here already, descendant has an index
										if(layerIndex[k]<=act_index){
											//we have been here already, suspicious for loops: TODO: this does not detect loops
											if(act_index>=modelbuf.size()){
												//the concurrency is endless: circular dependency
												printf("[Error]: There is circular dependency between models.\n");
												printf("         Models %s and %s are parts of a loop.\n", descendant->modelId().c_str(), modelbuf[i]->modelId().c_str());
												mTreeError=true;
												return;
											}
											else{
												layerIndex[k]=act_index+1;	//concurrent paths, apply longest so far
												if(max_layer_index<act_index+1){
													max_layer_index=act_index+1;
												}
											}
										}
									}
									else{
										numberofdescendants++;
										layerIndex[k]=act_index+1;
										break;
									}
								}
							}
						}
					}
				}
			}
			if(max_layer_index<act_index){
				max_layer_index=act_index;
			}
			act_index++;
		}while(numberofdescendants>0);
	}
	
	//warn if there are still unclassified models
	bool firstwarn=true;
	for(int j=0; j<layerIndex.size(); j++){
		if(layerIndex[j]==-1){
			if(firstwarn){
				printf("[Error]: There are models outside the network hierarchy:\n");
				firstwarn=false;
			}
			printf("\t#%d\t%s\n",j,mModels[j]->modelId().c_str());
		}
	}
	
	mModels.clear();
	
	for(int i=max_layer_index; i>=0; i--){
		for(int j=0; j<modelbuf.size(); j++){
			if(layerIndex[j]==i){
				mModels.push_back(modelbuf[j]);	//need to make them unconst for solving
			}
		}
	}
	
	//DEBUG
	//printf("\n*** SOLUTION ORDER ***\n");
	for(int i=0; i<mModels.size(); i++){
		//printf("\t%d\t%s\n",i+1,mModels[i]->modelId().c_str());
	}
	//printf("\n");
	//END DEBUG
}

//--------------------------------------------------------------------------------------------------

bool iWQSolver::saveInitVals(iWQInitialValues * yfrom)
{			
	if(!yfrom){
		return false;
	}
	int i=0;
	for(i=0; i<mModels.size(); i++){
        //for each valid model set the initial values
        if(mModels[i]){
        	mModels[i]->setInitialValues(yfrom);
        }
    }
    for(i=0; i<mExportLinks.size(); i++){
		mExportLinks[i].zerodest();
	}
	for(i=0; i<mExportLinks.size(); i++){
		mExportLinks[i].linkadd();
	}
	return true;
}

//--------------------------------------------------------------------------------------------------

bool iWQSolver::solve1Step(double xfrom, double xto, iWQInitialValues * yfrom)
{			
	if(xto<xfrom){
		//DEBUG
		//printf("Invalid solution order from %lf to %lf.\n",xfrom, xto);
		//END DEBUG
		return true;
	}
	
	int i, j;
	bool cleansolution=true;
	if(yfrom){	//forget wrong models only in the beginning
		mFaultyModels.clear();
	}
	
	for(i=0; i<mLinks.size(); i++){
		mLinks[i].zerodest();
	}
	for(i=0; i<mLinks.size(); i++){
		mLinks[i].linkadd();
	}
	
	for(i=0; i<mModels.size(); i++){
        //solve models in dependency order
		if(!mModels[i]->solve1Step(xfrom, xto, yfrom, mHmin, mEps)){
			mFaultyModels.push_back(mModels[i]);
			cleansolution=false;
		}
		for(j=0; j<mInterLinks.size(); j++){
			mInterLinks[j].zerodest();
		}
		for(j=0; j<mInterLinks.size(); j++){
			mInterLinks[j].linkadd();
		}
	}
	
	for(i=0; i<mExportLinks.size(); i++){
		mExportLinks[i].zerodest();
	}
	for(i=0; i<mExportLinks.size(); i++){
		mExportLinks[i].linkadd();
	}
	return cleansolution;
}

//--------------------------------------------------------------------------------------------------

bool iWQSolver::valid()
{
	return !mTreeError;
}

//--------------------------------------------------------------------------------------------------

void iWQSolver::setMinStepLength(double value)
{
	if(value>0.0){
		mHmin=value;
	}
	else{
		printf("[Warning]: Invalid minimal step length value (%lf) for solver.\n",value);
	}
}

//--------------------------------------------------------------------------------------------------

void iWQSolver::setAccuracy(double value)
{
	if(value>0.0){
		mEps=value;
	}
	else{
		printf("[Warning]: Invalid accuracy value (%lf) for solver.\n",value);
	}
}

//--------------------------------------------------------------------------------------------------

std::vector<std::string> iWQSolver::exportedDataHeaders(iWQDataTable * datatable)
{
	//returns data headers of the outputs
	std::vector<std::string> result;
	
	if(!datatable){
		return result;
	}
	
	std::vector<double *> ptrs;
	for(int i=0; i<mExportLinks.size(); i++){
		double * ptr=mExportLinks[i].destptr;
		//check if it is already there (multiple links can point to the same field
		bool alreadythere=false;
		for(int j=0; j<ptrs.size(); j++){
			if(ptrs[j]==ptr){
				alreadythere=true;
				break;
			}
		}
		if(!alreadythere){
			ptrs.push_back(ptr);
			//now check the name
			std::string name=datatable->columnForPort(ptr);
			result.push_back(name);
		}
	}
	
	return result;
}

//--------------------------------------------------------------------------------------------------

std::map<std::string, iWQKeyValues> iWQSolver::modelState()
{
	//save the state of all models in a map
	std::map<std::string, iWQKeyValues> result;
	//bool firstwithoutid=true;
	for(int i=0; i<mModels.size(); i++){
		iWQModel * act_model=mModels[i];
		if(act_model){
			std::string id=act_model->modelId();
			iWQStrings varnames=act_model->variableNames();
			iWQKeyValues state;
			//load state
			for(int j=0; j<varnames.size(); j++){
				std::string act_varname=varnames[j];
				const double * outlet=act_model->routlet(act_varname);
				state[act_varname]=*outlet;
			}
			
			if(id.size()==0){
				/*if(firstwithoutid){
					//save them as default
					result[""]=state;
				}
				else{
					//warn about that
					printf("[Error]: Model #%d is not the only one without id. Its state could not be saved.\n");
				}*/
			}
			else{
				//has id
				result[id]=state;
			}
		}
	}
	return result;
}

//--------------------------------------------------------------------------------------------------

void iWQSolver::setModelState(std::map<std::string, iWQKeyValues> state)
{
	for(int i=0; i<mModels.size(); i++){
		iWQModel * act_model=mModels[i];
		if(act_model){
			act_model->resetState();
			std::string id=act_model->modelId();
			iWQStrings varnames=act_model->variableNames();
			std::map<std::string, iWQKeyValues>::iterator it=state.find(id);
			if(it!=state.end()){
				iWQKeyValues act_state=it->second;
				iWQKeyValues::iterator it2;
				for(int j=0; j<varnames.size(); j++){
					it2=act_state.find(varnames[j]);
					double val=0.0;
					if(it2!=act_state.end()){
						//has value for that
						val=it2->second;
					}
					//set value 
					act_model->setStateVariable(varnames[j],val);
				}
			}
			else{
				printf("[Error]: No state information for %s.\n",id.c_str());
			}
		}
		else{
			printf("[Error]: Invalid model pointer.\n");
		}
	}
}

//--------------------------------------------------------------------------------------------------

std::vector<iWQModel *> iWQSolver::modelsThatDidNotSolve()
{
	return mFaultyModels;
}


