/*
 *  setup.cpp
 *  Highest-level object connecting everything and interpreting layouts
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/SETUP
 *
 */ 

#include <stdio.h>
#include <map>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <sys/stat.h>
#include <time.h>
#include <float.h>

#include "tinyxml/tinyxml.h"
#include "setup.h"
#include "model.h"
#include "solver.h"
#include "evaluator.h"
#include "datatable.h"
#include "modelfactory.h"
#include "mathutils.h"
#include "evaluatormethod.h"
#include "sampleutils.h"
#include "seriesinterface.h"
#include "filter.h"
#include "script.h"

//BEGIN NEW
#include "Eigen/Dense"
#include "biasmatrices.h"
//END NEW

#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

#define IWQ_LAYOUT_MIN_VERSION 0.2

//---------------------------------------------------------------------------------------

void SaveMatrix(Eigen::MatrixXd matrix, std::string filename);
Eigen::MatrixXd LoadMatrix(std::string filename);

//---------------------------------------------------------------------------------------

#ifdef _WIN32
//unix utils which are already in libmodel.a
extern int usleep(int useconds);
extern int mkdir(const char * name, int accessflags);
#endif

void iWQModelLayout::printError(std::string errormessage, TiXmlElement * element, int errorlevel)
{
	std::string fileposstr="";
	std::string errorlabel;
	if(errorlevel<=0){
		errorlabel="Warning";
	}
	else{
		errorlabel="Error";
	}
	if(element){
		std::stringstream s (fileposstr);
		s<<" "<<mFilename;
		s<<":";
		s<<element->Row();
		fileposstr=s.str();
	}
	printf("[%s]:%s %s\n",errorlabel.c_str(),fileposstr.c_str(),errormessage.c_str());
}

//---------------------------------------------------------------------------------------

iWQModelLayout::iWQModelLayout(std::string filename)
{
	mModels.clear();
	mLinks.clear();
	mExportLinks.clear();
	mSolver=NULL;
	mCommonParameters=NULL;
	mEvaluator=NULL;
	//mEvaluatorMethod=NULL;
	mDataTable=NULL;
	mComparisonLinks.clear();
	mInitVals=NULL;
	mFilename="";
	mSeriesInterface=NULL;
	mFilters.clear();
	mPreScripts.clear();
	mPostScripts.clear();
	
	//initialize model factory
	mModelFactory=new iWQModelFactory("models");
	
	//read conf from an xml file
	TiXmlDocument * doc=new TiXmlDocument(filename.c_str());
	if(!doc->LoadFile()){
		printf("[Error]: Failed to load XML model description file.\n");
		if(doc->Error()){
			printf("%s:%d (%d) error code %d:\t%s\n",filename.c_str(),doc->ErrorRow(),doc->ErrorCol(),doc->ErrorId(),doc->ErrorDesc());
		}
		return;
	}
	else{
		mFilename=filename;
	}
	if(doc->Error()){
		printf("[Error]: %s:%d %s.\n",filename.c_str(),doc->ErrorRow(), doc->ErrorDesc());
		return;
	}
	
	TiXmlHandle docHandle (doc);
	
	//check layout version
	if(!checkLayoutVersion(docHandle)){
		printError("[version] of <layout> not supported (should be above " QUOTEME(IWQ_LAYOUT_MIN_VERSION) ").",NULL);
	}
	else{
		//MODELS
		loadModels(docHandle);
		
		//DATA (only the 1st node of this type is processed)
		loadData(docHandle);
		
		//DISTRIBUTIONS
		loadDistributions(docHandle);
		
		//PARAMETERS
		loadParameters(docHandle);
		
		//CONNECTIONS
		loadConnections(docHandle);
		
		//INITIAL VALUES
		loadInitVals(docHandle);
		
		//COMPARISON LINKS & EVALUATOR
		loadComparisonLinks(docHandle);

		//FILTERS
		loadFilters(docHandle);
		
		//SCRIPTS
		loadScripts(docHandle);
		
		//INIT SOLVER
		//if(mExportLinks.size()){
			mSolver=new iWQSolver(mLinks,mExportLinks);
		//}
		
		if(!mSolver || !mSolver->valid()){
			printf("[Error]: Could not create solver.\n");
		}
		
		//CONFIG SOLVER
		configureSolver(docHandle);
		
		//LOAD EVALUATION METHOD
		//loadEvaluationMethod(docHandle); //anyway called by loadComparisonLinks
		
		//INIT EVALUATOR
		if(mDataTable && mSolver && mSolver->valid() && mCommonParameters && mInitVals && mComparisonLinks.size()>0){
			mEvaluator=new iWQEvaluator;
			mEvaluator->setDataTable(mDataTable);
			mEvaluator->setSolver(mSolver);
			mEvaluator->setParameters(mCommonParameters);
			mEvaluator->setInitialValues(mInitVals);
			mEvaluator->setComparisonLinks(mComparisonLinks);
			mEvaluator->setEvaluatorMethods(mEvaluatorMethods);
			mEvaluator->setEvaluatorWeights(mEvaluatorWeights);
			mEvaluator->setFilters(mFilters);
			mEvaluator->setPreScripts(mPreScripts);
			mEvaluator->setPostScripts(mPostScripts);
		}
		
		//CONFIG OPTIMIZER
		configureOptimizer(docHandle);
	}
	
	delete doc;
}

//---------------------------------------------------------------------------------------

iWQModelLayout::~iWQModelLayout()
{
	if(mEvaluator) delete mEvaluator;
	if(mSeriesInterface) delete mSeriesInterface;
	if(mDataTable) delete mDataTable;
	if(mSolver) delete mSolver;
	if(mCommonParameters) delete mCommonParameters;
	if(mModelFactory){
		for(unsigned int i=0; i<mModels.size(); i++){
			if(mModels[i]) mModelFactory->deleteModel(mModels[i]); //delete mModels[i];
		}
		delete mModelFactory;
	}
	//delete distributions
	std::map<std::string, iWQDistribution * >::iterator it;
	for(it=mDistributions.begin(); it!=mDistributions.end(); ++it){
		delete it->second;
	}
}

//---------------------------------------------------------------------------------------

#pragma mark XML parsing

void iWQModelLayout::loadModels(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xmodel=docHandle.FirstChild("layout").FirstChild("model").ToElement();
	while(xmodel){
		//process it
		std::string modeltype="";
		std::string modelid="";
		std::vector<std::string> modelflags;
		std::map<std::string, double> ownparams;
		
		bool valid=true;
		if(xmodel->QueryStringAttribute("type",&modeltype)!=TIXML_SUCCESS){
			printError("<model> does not have a [type] attribute.",xmodel);
			valid=false;
		}
		if(xmodel->QueryStringAttribute("id",&modelid)!=TIXML_SUCCESS){
			printError("<model> does not have an [id] attribute.\n",xmodel,0);
		}
		//get attributes
		TiXmlNode * xmodelattribnode=xmodel->FirstChild("attribute");
		if(xmodelattribnode){
			TiXmlElement * xmodelattrib=xmodelattribnode->ToElement();
			while(xmodelattrib){
				std::string attribval=xmodelattrib->GetText();
				modelflags.push_back(attribval);
			
				//jump tot next
				next=xmodelattrib->NextSibling("attribute");
				if(next){
					xmodelattrib=next->ToElement();
				}
				else{
					break;
				}
			}
		}
		//get own parameters
		TiXmlNode * xmodelparamnode=xmodel->FirstChild("parameter");
		if(xmodelparamnode){
			TiXmlElement * xmodelparam=xmodelparamnode->ToElement();
			while(xmodelparam){
				std::string paramname="";
				double paramval=0.0;
				bool valid=true;
				if(xmodelparam->QueryStringAttribute("name",&paramname)!=TIXML_SUCCESS){
					valid=false;
				}
				if(xmodelparam->QueryDoubleAttribute("value", &paramval)){
					valid=false;
				}
				
				if(valid){
					ownparams[paramname]=paramval;
				}
								
				//jump tot next
				next=xmodelparam->NextSibling("parameter");
				if(next){
					xmodelparam=next->ToElement();
				}
				else{
					break;
				}
			}
		}
		//we have all the values
		if(valid){
			
			iWQModel * model=mModelFactory->newModelOfType(modeltype);;
			
			//assemble model
			if(model){
				model->setModelId(modelid);
				model->setModelFlags(modelflags); 
				std::map<std::string, double>::iterator it;
				for(it=ownparams.begin(); it!=ownparams.end(); ++it){
					model->setValueForParam(it->second, it->first);
				}
				mModels.push_back(model);
			}
		}
		else{
			printError("Invalid <model> node.",xmodel);
		}
		
		//jump to next
		next=xmodel->NextSibling("model");
		if(next){
			xmodel=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadParameters(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xpar=docHandle.FirstChild("layout").FirstChild("parameters").ToElement();
	if(xpar){
		std::string parfilename="";
		mCommonParameters=new iWQParameterManager;
		if(xpar->QueryStringAttribute("src",&parfilename)==TIXML_SUCCESS){
			if(parfilename.size()){
				mCommonParameters->initFromFile(parfilename);
			}
			else{
				printError("Empty [src] attribute specified in <parameters>.",xpar);
			}
			//warn if there are also explicit values
			if(!xpar->NoChildren()){
				printError("Explicit <parameter> declarations are skipped if an [src] attribute is present in <parameters>.",xpar,0);
			}
		}
		else{
			//look foor explicit declarations
			TiXmlNode * xexpparamnode=xpar->FirstChild("parameter");
			if(xexpparamnode){
				TiXmlElement * xexpparam=xexpparamnode->ToElement();
				while(xexpparam){
					std::string paramname="";
					std::string domain="";
					double paramval=0.0;
					std::string distribution="";	//prior dist name 
					bool valid=true;
					if(xexpparam->QueryStringAttribute("name",&paramname)!=TIXML_SUCCESS){
						valid=false;
					}
					if(xexpparam->QueryDoubleAttribute("value", &paramval)!=TIXML_SUCCESS){
						valid=false;
					}
					xexpparam->QueryStringAttribute("domain", &domain);
					xexpparam->QueryStringAttribute("distribution", &distribution);
				
					if(valid){
						mCommonParameters->initParam(paramname, domain, paramval);
						//try to attach the given distribution
						if(distribution.size()){
							std::map<std::string, iWQDistribution * >::iterator it=mDistributions.find(distribution);
							if(it!=mDistributions.end()){
								mCommonParameters->linkDistributionToParam(it->second, paramname, domain);
							}
							else{
								printError("Unknown [distribution] name.",xexpparam);
							}
						}
					}
					else{
						printError("Invalid explicit <parameter> declaration.",xexpparam);
					}
					
					//check limits
					iWQLimits lim;
					bool hasmin=false;
					bool hasmax=false;
					if(xexpparam->QueryDoubleAttribute("min", &(lim.min))==TIXML_SUCCESS){
						hasmin=true;
					}
					if(xexpparam->QueryDoubleAttribute("max", &(lim.max))==TIXML_SUCCESS){
						hasmax=true;
					}
								
					//set limits (if applicable)
					if(hasmin && hasmax && (paramval>=lim.min && paramval<=lim.max) && lim.min<=lim.max){
						mCommonParameters->setLimitsForParam(lim, paramname, domain);
					}
					else{
						//wrong specification
						if(hasmin || hasmax){
							printError("<parameter> limits should contain both [min] and [max].",xexpparam);
						}
						else{
							if(hasmin && hasmax){
								if(paramval<lim.min || paramval>lim.max){
									printError("The [value] of <parameter> is outside the limits.",xexpparam);
								}
								if(lim.min>lim.max){
									printError("[min] should be less or equal to [max] in <parameter>.",xexpparam);
								}
							}
						}
					}
					
					//jump tot next
					next=xexpparam->NextSibling("parameter");
					if(next){
						xexpparam=next->ToElement();
					}
					else{
						break;
					}
				}
			}
		}
		for(unsigned int i=0; i<mModels.size(); i++){
			mModels[i]->bind(mCommonParameters);
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadData(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xdata=docHandle.FirstChild("layout").FirstChild("data").ToElement();
	
	std::map<std::string,std::string> inputserieslinks;	//these hold the linked fields temporarily
	std::map<std::string,std::string> outputserieslinks;
	
	if(xdata){
		std::string datafilename="";
		mDataColsToExport.clear();
		if(xdata->QueryStringAttribute("src",&datafilename)==TIXML_SUCCESS){
			mDataTable=new iWQDataTable (datafilename);
			//check if loaded correctly
			if(!mDataTable->numRows()){
				//empty or corrupted file
				printError("Data table could not be initialized from file \""+datafilename+"\"",xdata);
				delete mDataTable;
				mDataTable=NULL;
			}
			//process new columns
			TiXmlNode * xnewcolnode=xdata->FirstChild("column");
			if(xnewcolnode && mDataTable){
				TiXmlElement * xnewcol=xnewcolnode->ToElement();
				while(xnewcol){
					//create it
					std::string newcol=xnewcol->GetText();
					mDataTable->addColumn(newcol);
					
					//add export directive
					std::string exportitstr;
					if(xnewcol->QueryStringAttribute("export",&exportitstr)==TIXML_SUCCESS){
						//ignore case
						std::transform(exportitstr.begin(), exportitstr.end(), exportitstr.begin(), ::tolower);
						bool exportit=false;
						//printf("[Debug]: export=\"%s\" for %s\n",exportitstr.c_str(), newcol.c_str());
						if(exportitstr.compare("1")==0 || exportitstr.compare("true")==0){
							exportit=true;
						}
						if(exportit){
							mDataColsToExport.push_back(newcol);
						}
					}
					//read series links
					std::string seriesfilestr;
					bool imported=false;
					if(xnewcol->QueryStringAttribute("srcseries",&seriesfilestr)==TIXML_SUCCESS){
						inputserieslinks[newcol]=seriesfilestr;
						imported=true;
					}
					if(xnewcol->QueryStringAttribute("destseries",&seriesfilestr)==TIXML_SUCCESS){
						if(imported){
								printError("<column> must have either an [srcseries] or a [destseries] attribute, but not both.",xdata);
						}
						else{
							outputserieslinks[newcol]=seriesfilestr;
						}
					}
			
					//jump to next
					next=xnewcol->NextSibling("column");
					if(next){
						xnewcol=next->ToElement();
					}
					else{
						break;
					}
				}
			}
			//set time column
			std::string timecolname="";
			if(xdata->QueryStringAttribute("timecol",&timecolname)==TIXML_SUCCESS){
				mDataTable->setTField(timecolname);
				if(!mDataTable->timePort()){
					timecolname="";
				}
			}
			if(timecolname.size()==0){			
				mDataTable->setTField(0);
			}
			//Now load the series interface
			mSeriesInterface=new iWQSeriesInterface(mDataTable);
			std::map<std::string, std::string>::iterator it;
			for(it=inputserieslinks.begin(); it!=inputserieslinks.end(); ++it){
				mSeriesInterface->addSeriesLink(it->first, it->second, false);
			}
			for(it=outputserieslinks.begin(); it!=outputserieslinks.end(); ++it){
				mSeriesInterface->addSeriesLink(it->first, it->second, true);
			}
		}
		else{
			printError("<data> must have an [src] attribute.",xdata);
		}
		if(xdata->NextSibling("data")){
			printError("Only the first <data> tag is processed.",xdata,0);
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadConnections(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xconn=docHandle.FirstChild("layout").FirstChild("connection").ToElement();
	while(xconn){
		bool fromdata=false;
		bool frommodel=false;
		bool todata=false;
		bool tomodel=false;
		bool totype=false;
		std::string sourcedata="";
		std::string destdata="";
		std::string sourceobj="";
		std::string sourceport="";
		std::string destobj="";
		std::string destport="";
		std::string desttype="";
		bool fixprop=false;
		bool keyprop=false;
		double proportion=1.0;
		std::string propkey="";
		
		//get attributes
		if(xconn->QueryStringAttribute("sourcedata",&sourcedata)==TIXML_SUCCESS){
			fromdata=true;
		}
		if(xconn->QueryStringAttribute("destdata",&destdata)==TIXML_SUCCESS){
			todata=true;
		}
		if(xconn->QueryStringAttribute("sourceobj",&sourceobj)==TIXML_SUCCESS && xconn->QueryStringAttribute("sourceport",&sourceport)==TIXML_SUCCESS){
			frommodel=true;
		}
		if(xconn->QueryStringAttribute("destobj",&destobj)==TIXML_SUCCESS && xconn->QueryStringAttribute("destport",&destport)==TIXML_SUCCESS){
			tomodel=true;
		}
		if(xconn->QueryStringAttribute("desttype",&desttype)==TIXML_SUCCESS && xconn->QueryStringAttribute("destport",&destport)==TIXML_SUCCESS){
			totype=true;
		}
		if(xconn->QueryDoubleAttribute("proportion", &proportion)==TIXML_SUCCESS){
			fixprop=true;
		}
		if(xconn->QueryStringAttribute("proportionkey", &propkey)==TIXML_SUCCESS){
			keyprop=true;
		}
		
		//filter for errors
		if(frommodel && fromdata){
			printError("<connection> should have either [sourcedata] or [sourceobj]+[sourceport] attributes, but never both.",xconn);
		}
		else if((tomodel && todata) || (tomodel && totype) || (todata && totype)){
			printError("<connection> should have either [destdata] or [destobj]+[destport] or [desttype]+[destport] attributes, but only one.",xconn);
		}
		else if(!frommodel && !fromdata){
			printError("<connection> should have [sourcedata] or [sourceobj]+[sourceport] attributes.",xconn);
		}
		else if(!tomodel && !todata && !totype){
			printError("<connection> should have [destdata] or [destobj]+[destport] or [desttype]+[destport] attributes.",xconn);
		}
		else{
			// seems to be valid connection (in terms of source and dest definitions)
			iWQLink link;
			double * psrcport=NULL;
			double * pdestport=NULL;
			iWQModel * psrcmodel=NULL;
			iWQModel * pdestmodel=NULL;
			std::vector<iWQModel *> pdestmodels;
			
			if(fromdata && mDataTable){
				psrcport=mDataTable->portForColumn(sourcedata);
			}
			if(todata && mDataTable){
				pdestport=mDataTable->portForColumn(destdata);
			}
			if(frommodel){
				//look for its ID in mModels
				for(unsigned int i=0; i<mModels.size(); i++){
					std::string act_id=mModels[i]->modelId();
					if(act_id.compare(sourceobj)==0){
						psrcmodel=mModels[i];
						break;
					}
				}
				//check model and port validity
				if(!psrcmodel){
					printError("<connection> refers to an invalid [sourceobj].",xconn);
				}
				else{
					if(!psrcmodel->routlet(sourceport)){
						printError("<connection> refers to an invalid [sourceport].",xconn);
						psrcmodel=NULL;
					}
				}
				
			}
			if(tomodel){
				//look for its ID in mModels
				for(unsigned int i=0; i<mModels.size(); i++){
					std::string act_id=mModels[i]->modelId();
					if(act_id.compare(destobj)==0){
						pdestmodel=mModels[i];
						break;
					}
				}
				//check port validity
				if(!pdestmodel ){
					printError("<connection> points to an invalid [destobj].",xconn);
				}
				else{
					if(!pdestmodel->rwoutlet(destport)){
						printError("<connection> points to an invalid [destport].",xconn);
						pdestmodel=NULL;
					}
				}
			}
			if(totype){
				//look for types of desttype in mModels
				bool error=false;
				for(unsigned int i=0; i<mModels.size(); i++){
					iWQModel * act_model=mModels[i];
					std::string act_type=act_model->modelType();
					if(act_type.compare(desttype)==0){
						//matching object
						if(act_model->rwoutlet(destport)){
							pdestmodels.push_back(act_model);
						}
						else{
							error=true;
						}
					}
				}
				if(error){
					printError("<connection> points to an invalid [destport] in models of [desttype].",xconn);
				}
			}
			//check proportional definitions
			if(fixprop && keyprop){
				printError("<connection> has both [proportion] and [proportionkey], both omitted.",xconn);
				fixprop=false;
				keyprop=false;
			}
			if(keyprop && (fromdata || todata)){
				printError("[proportionkey] cannot be used in a <connection> attached to data, omitted.",xconn);
				fixprop=false;
				keyprop=false;
			}
			if(fixprop && proportion==0){
				printError("[proportion] is 0.",xconn,0);	//was: "[proportion] is 0, omitted".
				//fixprop=false;	//allow 0 proportions for switching off connections parametrically
				//keyprop=false;
			}
			if(fixprop && proportion<0){
				printError("<connection> has negative fixed [proportion].",xconn,0);
			}
			if(keyprop && propkey.size()==0){
				printError("[proportionkey] is empty, omitted.",xconn);
				fixprop=false;
				keyprop=false;
			}
						
			// pure data links
			if(fromdata && todata){
				printError("<connection> simply copies data.",xconn,0);
				if(!psrcport || !pdestport){
					printError("Invalid <connection> of type \"data->data\".",xconn);
				}
				else{
					link.establish(psrcport,pdestport);
					if(fixprop){
						link.setFixedProportion(proportion);
					}
					//load if not already there
					if(std::find(mExportLinks.begin(), mExportLinks.end(), link)==mExportLinks.end()){
						mExportLinks.push_back(link);
					}
					else{
						printError("This <connection> has been already defined elsewhere.",xconn,0);
					}
				}
			}
			//there is exactly 1 source and 1 destination
			if(fromdata && tomodel){
				if(psrcport && pdestmodel){
					link.establish(psrcport, pdestmodel, destport);
					if(fixprop){
						link.setFixedProportion(proportion);
					}
					//load if not already there
					if(std::find(mLinks.begin(), mLinks.end(), link)==mLinks.end()){
						mLinks.push_back(link);
					}
					else{
						printError("This <connection> has been already defined elsewhere.",xconn,0);
					}
				}
				else{
					printError("Invalid <connection> of type \"data->model\".",xconn);
				}				
			}
			if(frommodel && todata){
				if(psrcmodel && pdestport){
					link.establish(psrcmodel, sourceport, pdestport);
					if(fixprop){
						link.setFixedProportion(proportion);
					}
					//load if not already there
					if(std::find(mExportLinks.begin(), mExportLinks.end(), link)==mExportLinks.end()){
						mExportLinks.push_back(link);
					}
					else{
						printError("This <connection> has been already defined elsewhere.",xconn,0);
					}
				}
				else{
					printError("Invalid <connection> of type \"model->data\".",xconn);
				}				
			}
			if(frommodel && tomodel){
				if(psrcmodel && pdestmodel){
					link.establish(psrcmodel, sourceport, pdestmodel, destport);
					if(fixprop){
						link.setFixedProportion(proportion);
					}
					if(keyprop){
						const double * srcoutlet=psrcmodel->routlet(propkey);
						const double * destoutlet=pdestmodel->routlet(propkey);
						if(srcoutlet && destoutlet){
							link.setKeyedProportion(propkey);
						}
						else{
							printError("Dependent proportion key of <connection> is not valid in a connected model.",xconn);
						}
					}
					//load if not already there
					if(std::find(mLinks.begin(), mLinks.end(), link)==mLinks.end()){
						mLinks.push_back(link);
					}
					else{
						printError("This <connection> has been already defined elsewhere.",xconn,0);
					}
				}
				else{
					printError("Invalid <connection> of type \"model->model\".",xconn);
				}				
			}
			//global links
			if(fromdata && totype){
				if(psrcport && pdestmodels.size()){
					for(unsigned int i=0; i<pdestmodels.size(); i++){
						link.establish(psrcport, pdestmodels[i], destport);
						if(fixprop){
							link.setFixedProportion(proportion);
						}
						//load if not already there
						if(std::find(mLinks.begin(), mLinks.end(), link)==mLinks.end()){
							mLinks.push_back(link);	
						}
						else{
							printError("This <connection> has been already defined elsewhere.",xconn,0);
						}
					}
				}
				else{
					printError("Invalid <connection> of type \"data->model_type\".",xconn);
				}				
			}
			if(frommodel && totype){
				if(psrcmodel && pdestmodels.size()){
					bool errorshown=false;
					for(unsigned int i=0; i<pdestmodels.size(); i++){
						link.establish(psrcmodel, sourceport, pdestmodels[i], destport);
						if(fixprop){
							link.setFixedProportion(proportion);
						}
						if(keyprop){
							const double * srcoutlet=psrcmodel->routlet(propkey);
							const double * destoutlet=pdestmodels[i]->routlet(propkey);
							if(srcoutlet && destoutlet){
								link.setKeyedProportion(propkey);
							}
							else{
								if(!errorshown){
									printError("Dependent proportion key of <connection> is not valid in a connected model.",xconn);
									errorshown=true;
								}
							}
						}
						//load if not already there
						if(std::find(mLinks.begin(), mLinks.end(), link)==mLinks.end()){
							mLinks.push_back(link);	
						}
						else{
							printError("This <connection> has been already defined elsewhere.",xconn,0);
						}
					}
				}
				else{
					printError("Invalid <connection> of type \"model->model_type\".",xconn);
				}
			}
		}
		
		//jump to next
		next=xconn->NextSibling("connection");
		if(next){
			xconn=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadInitVals(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xinit=docHandle.FirstChild("layout").FirstChild("initials").FirstChild("initial").ToElement();
	mInitVals=new iWQInitialValues;
	while(xinit){
				
		//process it
		std::string varname="";
		std::string modelid="";
		double value=0.0;
			
		bool valid=true;
		if(xinit->QueryStringAttribute("variable",&varname)!=TIXML_SUCCESS){
			printError("<initial> does not have a [variable] attribute.",xinit);
			valid=false;
		}
		if(xinit->QueryDoubleAttribute("value",&value)!=TIXML_SUCCESS){
			printError("<initial> does not have a [value] attribute.",xinit);
			valid=false;
		}
		if(valid){
			xinit->QueryStringAttribute("model",&modelid);
			if(modelid.size()){
				mInitVals->setValueForVariable(value, varname, modelid);
			}
			else{
				mInitVals->setDefaultValueForVariable(value, varname);
			}
		}		
		
		//jump to next
		next=xinit->NextSibling("initial");
		if(next){
			xinit=next->ToElement();
		}
		else{
			break;
		}
	}
	
	//wire up with the parameters
	mInitVals->setParameterManager(mCommonParameters);
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadFilters(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xfilter=docHandle.FirstChild("layout").FirstChild("filters").FirstChild("filter").ToElement();
	mFilters.clear();
	while(xfilter){
				
		//process it
		std::string srcname="";
		std::string destname="";
		int winlen = -1;
		int wincenter = -1;
		std::string function = "";
			
		bool valid=true;	//xml validity
		if(xfilter->QueryStringAttribute("sourcedata",&srcname)!=TIXML_SUCCESS){
			printError("<filter> does not have a [sourcedata] attribute.",xfilter);
			valid=false;
		}
		if(xfilter->QueryStringAttribute("destdata",&destname)!=TIXML_SUCCESS){
			printError("<filter> does not have a [destdata] attribute.",xfilter);
			valid=false;
		}
		if(xfilter->QueryIntAttribute("length",&winlen)!=TIXML_SUCCESS){
			printError("<filter> does not have a [length] attribute.",xfilter);
			valid=false;
		}
		if(xfilter->QueryIntAttribute("center",&wincenter)!=TIXML_SUCCESS){
			printError("<filter> does not have a [center] attribute.",xfilter);
			valid=false;
		}
		if(xfilter->QueryStringAttribute("function",&function)!=TIXML_SUCCESS){
			printError("<filter> does not have a [function] attribute.",xfilter);
			valid=false;
		}
		if(valid){
			iWQFilter * f = new iWQFilter;
			if(f){
				if(	f->setDataTable(mDataTable) && 
					f->setSrcFieldName(srcname) && 
					f->setDestFieldName(destname) &&
					f->setFunction(function) &&
					f->setWindowLength(winlen) &&
					f->setWindowCenter(wincenter)){
					//valid data
					mFilters.push_back(f);
				}
			}
			
		}
				
		//jump to next
		next=xfilter->NextSibling("filter");
		if(next){
			xfilter=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadComparisonLinks(TiXmlHandle docHandle)
{
	//get the evaluation method name
	TiXmlNode * next;
	TiXmlElement * xevalmethod=docHandle.FirstChild("layout").FirstChild("evaluation").ToElement();
	while(xevalmethod){
		std::string evalmethodname="";
		if(xevalmethod->QueryStringAttribute("method",&evalmethodname)!=TIXML_SUCCESS){
			printError("<evaluation> should have a [method] attribute.",xevalmethod);
			return;
		}
		if(evalmethodname.size()==0){
			printError("The [method] attribute of <evaluation> should not be empty.",xevalmethod);
			return;
		}
		iWQEvaluatorMethod * testmethod=createEvalMethod(evalmethodname);
		if(!testmethod){
			printError("The [method] attribute of <evaluation> refers to an unknown error model.",xevalmethod);
			return;
		}
		else{
			delete testmethod;
		}
		
		//check if there are other <evaluation> tags
		//jump to next
		/*next=xevalmethod->NextSibling("evaluation");
		if(next){
			printError("Only the first <evaluation> tag is processed.",next->ToElement(),0);
		}*/
		
		TiXmlNode * xcomplinknode = xevalmethod->FirstChild("compare");
		if(!xcomplinknode){
			printError("<evaluation> should have at least one <compare> tag.",xevalmethod);
			return;
		}
		TiXmlElement * xcomplink = xcomplinknode->ToElement();
		while(xcomplink){
			//process it
			
			std::string modelledcol="";
			std::string measuredcol="";
			double weight=1.0;
			
			bool valid=true;
			if(xcomplink->QueryStringAttribute("modelled",&modelledcol)!=TIXML_SUCCESS){
				printError("<compare> should have a [modelled] attribute.",xcomplink);
				valid=false;
			}
			if(modelledcol.size()==0){
				printError("The [modelled] attribute of <compare> should not be empty.",xcomplink);
				valid=false;
			}
			if(xcomplink->QueryStringAttribute("measured",&measuredcol)!=TIXML_SUCCESS){
				printError("<compare> should have a [measured] attribute.",xcomplink);
				valid=false;
			}
			if(measuredcol.size()==0){
				printError("The [measured] attribute of <compare> should not be empty.",xcomplink);
				valid=false;
			}
			if(xcomplink->QueryDoubleAttribute("weight",&weight)==TIXML_SUCCESS){
				if(weight<=0.0){
					printError("This <compare> is omitted due invalid or zero [weight] attribute.",xcomplink,0);
					valid=false;
				}
			}
			if(valid){
				iWQComparisonLink cl (mDataTable, modelledcol, measuredcol);
				if(cl.valid()){
					//load if not already there
					if(std::find(mComparisonLinks.begin(),mComparisonLinks.end(),cl)==mComparisonLinks.end()){
						//load the corresponding evaluator method
						iWQEvaluatorMethod * method=loadEvaluationMethod(xcomplink,cl,evalmethodname);
						if(!method){
							printError("Failed to create evaluator method for <compare> tag.",xcomplink);
							return;
						}
						mEvaluatorMethods.push_back(method);
						mEvaluatorWeights.push_back(weight);
						//everything OK, load the comparison link too
						mComparisonLinks.push_back(cl);
					}
					else{
						printError("This <compare> tag has been already defined elsewhere.",xcomplink,0);
					}
				}
				else{
					printError("Invalid <compare> tag.",xcomplink);
				}
			}		
			
			//jump to next
			next=xcomplink->NextSibling("compare");
			if(next){
				xcomplink=next->ToElement();
			}
			else{
				break;
			}	
		}
	
		//jump to next evalmethod
		next=xevalmethod->NextSibling("evaluation");
		if(next){
			xevalmethod=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

iWQEvaluatorMethod * iWQModelLayout::loadEvaluationMethod(TiXmlElement * compareNode, iWQComparisonLink link, std::string methodName)
{
	if(compareNode==NULL){	// ||  || methodname.size()==0
		return NULL;
	}
	
	iWQEvaluatorMethod * evalMethod=createEvalMethod(methodName);
	std::multimap<std::string, std::string> othersettings;
	
	//now initialize evalMethod
	if(evalMethod){
		//set modelled col
		evalMethod->setComparisonLink(link);
		
		//read settings
		TiXmlNode * xsettings=compareNode->FirstChild("settings");
		bool hassettingsmap=(xsettings!=NULL);
		if(xsettings){
			//recursively store settings
			storeNodeInMap(xsettings, &othersettings, "", 0);
		}
			
		//load settings
		if(evalMethod && evalMethod->wantsParams()){
			evalMethod->setParameterStorage(mCommonParameters);	//needs to do this because of early parameter loading
			//in-place settings
			//if(hassettingsmap){
			evalMethod->setParams(othersettings);	//need to load an empty set too, because of the dynamic parameters
			//}
		}
	}
		
	return evalMethod;
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadDistributions(TiXmlHandle docHandle)
{
	//prior distribution storage
	TiXmlNode * next;
	TiXmlElement * xdistrib=docHandle.FirstChild("layout").FirstChild("distribution").ToElement();
	
	mDistributions.clear();
	
	while(xdistrib){
		
		std::string name="";
		std::string type="";
		std::map<std::string, double> settings;
		
		//load the name
		if(xdistrib->QueryStringAttribute("name",&name)!=TIXML_SUCCESS){
			printError("<distribution> should have a [name] attribute.",xdistrib);
		}
		
		//check for the uniqueness of name
		if(mDistributions.find(name)!=mDistributions.end()){
			printError("This [name] has been already used for another <distribution>.",xdistrib);
		}
		
		//get the type
		if(xdistrib->QueryStringAttribute("type",&type)!=TIXML_SUCCESS){
			printError("<distribution> should have a [type] attribute.",xdistrib);
		}
		
		//alloc the ditribution
		iWQDistribution * dist=NULL;
		if(type.compare("normal")==0){
			dist=new iWQRandomNormalGenerator();
		}
		else if(type.compare("lognormal")==0){
			dist=new iWQRandomLogNormalGenerator(1.0, 1.0); 	//fake parameters
		}
		else if(type.compare("t")==0){
			dist=new iWQRandomtGenerator();
		}
		else if(type.compare("uniform")==0){
			dist=new iWQRandomUniformGenerator();
		}
		else if(type.compare("exponential")==0){
			dist=new iWQRandomExpGenerator();
		}
		else if(type.compare("beta")==0){
			dist=new iWQRandomBetaGenerator();
		}
		else if(type.compare("sep")==0){
			dist=new iWQRandomSEPGenerator();
		}
		else if(type.compare("gamma")==0){
			dist=new iWQRandomGammaGenerator();
		}
		//todo: other types
		
		if(dist==NULL){
			printError("Unknown [type] specified for <distribution>.",xdistrib);
		}
		
		//now load the other settings from the attributes
		TiXmlAttribute* pAttrib=xdistrib->FirstAttribute();
		std::string keyname;
		double val;
		while (dist && pAttrib){
			keyname=pAttrib->Name();
			
			if(keyname.compare("name")!=0 && keyname.compare("type")!=0){
				if(pAttrib->QueryDoubleValue(&val)==TIXML_SUCCESS){
					settings[keyname]=val;
				}
				else{
					printError("<distribution> should only have number attributes besides [name] and [type]",xdistrib);
					printf    ("\t The value for [%s] is not a valid number.\n",keyname.c_str());
				}
			}
			pAttrib=pAttrib->Next();
		}

		if(dist){
			dist->initialize(settings);
		}
		
		//load the final distribution
		if(dist){
			mDistributions[name]=dist;
		}
		else{
			printError("<distribution> ignored due to errors.",xdistrib);
		}
				
		next=xdistrib->NextSibling("distribution");
		if(next){
			xdistrib=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadScripts(TiXmlHandle docHandle)
{
	//external scripts
	TiXmlNode * next;
	TiXmlElement * xscript=docHandle.FirstChild("layout").FirstChild("script").ToElement();
	
	mPreScripts.clear();
	mPostScripts.clear();
	
	//<script phase="PRE" order="4" command="" inputtable="" outputtable="" inputparams="" tabdelimitedparameters="" /> 
	while(xscript){
		//...
		std::string command="";
		std::string phasestr="";
		std::string intablename="";
		std::string outtablename="";
		std::string inparname="";
		bool tabdelim=false;
		int order = 99;
		int iphase = -1;	//-1: error, 0: pre, 1: post
		bool error=false;
		
		//load the command
		if(xscript->QueryStringAttribute("command",&command)!=TIXML_SUCCESS){
			printError("<script> must have a [command] attribute.",xscript);
			error=true;
		}
		
		//load the phase
		if(xscript->QueryStringAttribute("phase",&phasestr)!=TIXML_SUCCESS){
			printError("<script> must have a [phase] attribute.",xscript);
			error=true;
		}
		else{
			std::transform(phasestr.begin(), phasestr.end(), phasestr.begin(), ::tolower);
			//try to decode the phase attribute
			if(phasestr.compare("pre")==0){
				iphase=0;
			}
			if(phasestr.compare("post")==0){
				iphase=1;
			}
		}
		
		//load the intablename
		if(xscript->QueryStringAttribute("inputtable",&intablename)!=TIXML_SUCCESS){
			printError("<script> should have an [inputtable] attribute.",xscript);
			error=true;
		}
		
		//load the outtablename
		if(xscript->QueryStringAttribute("outputtable",&outtablename)!=TIXML_SUCCESS){
			printError("<script> should have an [outputtable] attribute.",xscript);
			error=true;
		}
		
		//load the inparname
		if(xscript->QueryStringAttribute("inputparams",&inparname)!=TIXML_SUCCESS){
			printError("<script> should have an [inputparams] attribute.",xscript);
			error=true;
		}
		
		//load the order
		if(xscript->QueryIntAttribute("order",&order)!=TIXML_SUCCESS){
			printError("<script> should have an [order] attribute.",xscript);
			error=true;
		}
		
		std::string tabdelimstr;
		if(xscript->QueryStringAttribute("tabdelimitedparameters",&tabdelimstr)==TIXML_SUCCESS){
			std::transform(tabdelimstr.begin(), tabdelimstr.end(), tabdelimstr.begin(), ::tolower);
			if(tabdelimstr.compare("1")==0 || tabdelimstr.compare("true")==0){
				tabdelim=true;
			}
		}
		
		if(!error){
			//alloc the script
			iWQScript s;
			s.setCommandString(command);
			s.setExportTableName(intablename);
			s.setImportTableName(outtablename);
			s.setExportParametersName(inparname);
			s.setOrder(order);
			s.setDataTable(mDataTable);
			s.setParameterManager(mCommonParameters);
			s.setExportTabDelimitedParameters(tabdelim);
			if(iphase==0){
				mPreScripts.push_back(s);
			}
			if(iphase==1){
				mPostScripts.push_back(s);
			}
		}
		else{
			printError("<script> ignored due to errors.",xscript);
		}
				
		next=xscript->NextSibling("script");
		if(next){
			xscript=next->ToElement();
		}
		else{
			break;
		}
	}
	
	//sort scripts
	if(mPreScripts.size()){
		std::sort(mPreScripts.begin(), mPreScripts.end());
	}
	if(mPostScripts.size()){
		std::sort(mPostScripts.begin(), mPostScripts.end());
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::storeNodeInMap(TiXmlNode * pParent, std::multimap<std::string, std::string> * container, std::string prefix, int level)
{
	if(!pParent || !container){
		return;
	}
		
	TiXmlNode * pChild;
	TiXmlText * pText;
	int t = pParent->Type();
	int num;

	if(level>0){
	if(prefix.size()){
		prefix+=":";
	}

	if(t==TiXmlNode::TINYXML_ELEMENT){
		TiXmlElement * xelement=pParent->ToElement();
		prefix+=pParent->ValueStr();
		const char * txt=xelement->GetText();
		std::string value=(txt!=NULL)?txt:"";
		
		//attributes
		TiXmlAttribute* pAttrib=xelement->FirstAttribute();
		int attribindex=0;
		std::string rootprefix=prefix;
		std::string attribstr="";
		while(pAttrib){
			std::string aprefix=pAttrib->Name();
			std::string avalue=pAttrib->ValueStr();
			if(aprefix.size() && avalue.size()){
				if(attribindex>0){
					attribstr+=",";
				}
				attribstr+=aprefix+"="+avalue;
				attribindex++;
			}
			pAttrib=pAttrib->Next();
		}
		if(attribindex>0){
			container->insert(std::pair<std::string,std::string>(prefix,attribstr));
			prefix+="("+attribstr+")";
		}
		
		if(value.size()){
			container->insert(std::pair<std::string,std::string>(prefix,value));
		}
	}
	}
	if(level<=128){
	for(pChild=pParent->FirstChild(); pChild!=0; pChild=pChild->NextSibling()){
		storeNodeInMap(pChild, container, prefix, level+1);
	}
	}
	else{
		printError("XML hierarchy levels are limited to 128 in <settings> - child nodes of "+pParent->ValueStr()+" skipped.\n",pParent->ToElement());
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::configureSolver(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xsolver=docHandle.FirstChild("layout").FirstChild("solver").ToElement();
	bool accset=false;
	bool stepset=false;
	
	if(xsolver && !mSolver){
		printError("There is no solver to configure.",xsolver);
		return;
	}
	
	while(xsolver){
		double precision;
		double minstep;
		if(xsolver->QueryDoubleAttribute("accuracy",&precision)==TIXML_SUCCESS){
			if(accset){
				printError("[accuracy] was already set for <solver>, now overriding.",xsolver,0);
			}
			mSolver->setAccuracy(precision);
			accset=true;
		}
		if(xsolver->QueryDoubleAttribute("minsteplength",&minstep)==TIXML_SUCCESS){
			if(stepset){
				printError("[minsteplength] was already set for <solver>, now overriding.",xsolver,0);
			}
			mSolver->setMinStepLength(minstep);
			stepset=true;
		}
		
		//jump to next
		next=xsolver->NextSibling("solver");
		if(next){
			xsolver=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::configureOptimizer(TiXmlHandle docHandle)
{
	TiXmlNode * next;
	TiXmlElement * xopt=docHandle.FirstChild("layout").FirstChild("optimizer").ToElement();
		
	if(xopt && !mEvaluator){
		printError("There is no optimizer to configure.",xopt);
		return;
	}
	
	while(xopt){
		//check PSO
		TiXmlNode * pso=xopt->FirstChild("particle-swarm");
		TiXmlElement * xpso=NULL;
		if(pso){
			xpso=pso->ToElement();
		}
		if(xpso){
			bool active=false;
			int numrounds;
			int swarmsize;
			int numidlerounds;
			std::string activestr;
			if(xpso->QueryStringAttribute("active",&activestr)==TIXML_SUCCESS){
				std::transform(activestr.begin(), activestr.end(), activestr.begin(), ::tolower);
				if(activestr.compare("1")==0 || activestr.compare("true")==0){
					active=true;
				}
				mEvaluator->PSOActive=active;
			}
			if(xpso->QueryIntAttribute("maxnumrounds",&numrounds)==TIXML_SUCCESS){
				mEvaluator->PSOMaxNumRounds=numrounds;
			}
			if(xpso->QueryIntAttribute("idlerounds",&numidlerounds)==TIXML_SUCCESS){
				mEvaluator->PSOMaxIdleRounds=numidlerounds;
			}
			if(xpso->QueryIntAttribute("size",&swarmsize)==TIXML_SUCCESS){
				mEvaluator->PSOSwarmSize=swarmsize;
			}
			if(active){
				printf("[optimizer]: Particle Swarm optimization is active (size: %d, rounds: %d, idlelimit: %d)\n",swarmsize,numrounds,numidlerounds);
			}
		}
		
		//check NMS
		TiXmlNode * nms=xopt->FirstChild("nelder-mead");
		TiXmlElement * xnms=NULL;
		if(nms){
			xnms=nms->ToElement();
		}
		if(xnms){
			bool active=true;
			int numrounds;
			double tolerance;
			std::string activestr;
			if(xnms->QueryStringAttribute("active",&activestr)==TIXML_SUCCESS){
				std::transform(activestr.begin(), activestr.end(), activestr.begin(), ::tolower);
				if(activestr.compare("1")==0 || activestr.compare("true")==0){
					active=true;
				}
				mEvaluator->NMSActive=active;
			}
			if(xnms->QueryIntAttribute("maxnumrounds",&numrounds)==TIXML_SUCCESS){
				mEvaluator->NMSMaxNumRounds=numrounds;
			}
			if(xnms->QueryDoubleAttribute("tolerance",&tolerance)==TIXML_SUCCESS){
				mEvaluator->NMSTolerance=tolerance;
			}
			if(active){
				printf("[optimizer]: Nelder-Mead Simplex optimization is active (max. rounds: %d, tolerance: %g)\n",numrounds,tolerance);
			}
		}
			
		//jump to next
		next=xopt->NextSibling("optimizer");
		if(next){
			xopt=next->ToElement();
		}
		else{
			break;
		}
	}
}

//---------------------------------------------------------------------------------------
bool iWQModelLayout::checkLayoutVersion(TiXmlHandle docHandle)
{
	TiXmlElement * xroot=docHandle.FirstChild("layout").ToElement();
	if(xroot){
		double version=-1.0;
		if(xroot->QueryDoubleAttribute("version",&version)==TIXML_SUCCESS){
			//printf("Version %lf detected.\n",version);
			return (version>=IWQ_LAYOUT_MIN_VERSION);
		}
	}
	return false;
}
//---------------------------------------------------------------------------------------

#pragma mark Quality (settings and model behaviour)

iWQModelLayoutValidity iWQModelLayout::validity()
{
	iWQModelLayoutValidity result=IWQ_NOT_VALID;
	
	//check the existence of components
	if(mDataTable && mSolver && mSolver->valid() && mCommonParameters && mInitVals && mDataTable->timePort() && mDataTable->numRows()){
		//criteria for running
		result=IWQ_VALID_FOR_RUN;
		if(mEvaluator && mComparisonLinks.size()){
			//criteria for evaluation
			result=IWQ_VALID_FOR_CALIBRATE;
		}
	}
	
	//TODO: test connections with a sample graph tour 
	
	return result;
}

//---------------------------------------------------------------------------------------

bool iWQModelLayout::verify()
{
	printf("Diagnosing models...\n");
	//check models 1 by 1
	for(unsigned int i=0; i<mModels.size(); i++){
		if(!(mModels[i]->verify())){
			printf("Verification failed.\n");
			return false;
		}
	}	
	printf("Passed\n");
	return true;
}

//---------------------------------------------------------------------------------------

#pragma mark Run, evaluate & calibrate

//private method without validity check: returns if the solution is stable
bool iWQModelLayout::runmodel(int * firsterrorrow, double * firsterrort)
{
	mDataTable->rewind();
	double * t=mDataTable->timePort();
	double prev_t = *t;
	iWQInitialValues * yfeed=mInitVals;
	bool stable=true;
	
	//run PRE scripts
	bool scriptsok=true;
	for(int s=0; s<mPreScripts.size(); s++){
		bool thisok = mPreScripts[s].execute();
		if(!thisok){
			printf("[Error]: Script \"%s\" failed to run correctly (return code=%d).\n",mPreScripts[s].commandString().c_str(),mPreScripts[s].returnStatus());
			scriptsok=false;
		}
	}
	
	//feed initvals into the data table
	mDataTable->rewind(); //not necessary here, just for making clear what happens in the next step
	mSolver->saveInitVals(yfeed);
	
	//run models
	while(mDataTable->stepRow()!=-1){
		
		if(!mSolver->solve1Step(prev_t, *t, yfeed)){
			if(stable && firsterrorrow){
				*firsterrorrow = mDataTable->pos();
			}
			if(stable && firsterrort){
				*firsterrort = prev_t;
			}
			stable=false;	
		}
				
		prev_t = *t;
		yfeed=NULL;
	}
	
	//run POST scripts
	for(int s=0; s<mPostScripts.size(); s++){
		bool thisok = mPostScripts[s].execute();
		if(!thisok){
			printf("[Error]: Script \"%s\" failed to run correctly (return code=%d).\n",mPostScripts[s].commandString().c_str(),mPostScripts[s].returnStatus());
			scriptsok=false;
		}
	}
	
	//run filters
	for(int i=0; i<mFilters.size(); i++){
		iWQFilter * f=mFilters[i];
		if(f){
			f->filter();
		}
	}
	
	return stable && scriptsok;
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::run()
{
	if(validity()<IWQ_VALID_FOR_RUN){
		printf("[Error]: Model layout is not suitable to run.\n");
		return;
	}
	if(!verify()){
		printf("[Error]: Model layout contains defects.\n");
		return;
	}
	int firsterrorrow = -1;
	double firsterrort = -DBL_MAX;
	if(!runmodel(&firsterrorrow, &firsterrort)){
		printf("[Warning]: Numerical stability could not be achieved with the minimal stepsize of %e.\n",mSolver->minStepLength());
		std::vector<std::string> parnames=mCommonParameters->namesForPlainValues();
		std::vector<double> parvalues=mCommonParameters->plainValues();
		for(int i=0; i<parnames.size(); i++){
			printf("%s=%g  ",parnames[i].c_str(), parvalues[i]);
		}
		printf("\n");
		std::vector<iWQModel *> wrongs = mSolver->modelsThatDidNotSolve();
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
}

//---------------------------------------------------------------------------------------

double iWQModelLayout::evaluate()
{
	if(validity()<IWQ_VALID_FOR_CALIBRATE){
		printf("[Error]: Model layout is not suitable to evaluate.\n");
		return 0.0;
	}
	if(!verify()){
		printf("[Error]: Model layout contains defects.\n");
		return 0.0;
	}
	return mEvaluator->evaluate();
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::calibrate()
{
	if(validity()<IWQ_VALID_FOR_CALIBRATE){
		printf("[Error]: Model layout is not suitable for calibration.\n");
		return;
	}
	if(!verify()){
		printf("[Error]: Model layout contains defects.\n");
		return;
	}
	mEvaluator->printWarnings=false;
	mEvaluator->calibrate();
	mEvaluator->printWarnings=true;
}

//---------------------------------------------------------------------------------------

#pragma mark File I/O wrappers 

void iWQModelLayout::saveParameters(std::string filename, bool tabdelimited)
{
	if(mCommonParameters){
		mCommonParameters->saveToFile(filename, tabdelimited);
	}
	else{
		printf("[Warning]: No parameters to save.\n");
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::loadParameters(std::string filename, bool tabdelimited)
{
	if(mCommonParameters){
		if(tabdelimited){
			mCommonParameters->initFromTabDelimitedFile(filename);
		}
		else{	
			mCommonParameters->initFromFile(filename);
		}
	}
	else{
		printf("[Warning]: No parameters to save.\n");
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::saveResults(std::string filename)
{
	//TODO: gather result/reference columns from iwqcomparisonlinks + time field from datatable
	if(mDataTable){
		mDataTable->writeToFile(filename);
	}
	else{
		printf("[Warning]: No data to save.\n");
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::saveResultsUNCSIM(std::string filename)
{
	//TODO: gather result/reference columns from iwqcomparisonlinks + time field from datatable
	if(mDataTable){
		if(mDataColsToExport.size()){
			mDataTable->saveUNCSIMFormatToFile(filename,mDataColsToExport);
		}
		else{
			mDataTable->saveUNCSIMFormatToFile(filename);
		}
	}
	else{
		printf("[Warning]: No data to save.\n");
	}
}

//---------------------------------------------------------------------------------------

#pragma mark DOT output

void iWQModelLayout::saveLayoutGraph(std::string dotfilename)
{
	//save network topology in .DOT file format
	FILE * ofile=fopen(dotfilename.c_str(),"w");
	if(!ofile){
		printf("[Error]: Could not open DOT file (\"%s\") for writing.\n", dotfilename.c_str());
		return;
	}
	fprintf(ofile,"digraph layout {\n");
	
	for(unsigned int i=0; i<mLinks.size(); i++){
		iWQModel * src=mLinks[i].dependsOn();
		iWQModel * dest=mLinks[i].subject();
		//skip data input links
		if(src && dest){
			std::string srcname=src->modelId();
			std::string destname=dest->modelId();
			fprintf(ofile,"\t%s -> %s;\n", srcname.c_str(), destname.c_str());
		}
	}
	for(unsigned int i=0; i<mExportLinks.size(); i++){
		double * port=mExportLinks[i].destinationPort();		// this is a data slot for a modelled field
		//indicate outputs considered in calibration
		iWQModel * src=mExportLinks[i].dependsOn();				// this is the modelled firld itself
		if(src){
			std::string srcname=src->modelId();
			std::string outputname=mDataTable->columnForPort(port);
			bool important=false;
			std::string meas_pair="";
			for(unsigned int j=0; j<mComparisonLinks.size(); j++){
				//see if that link applies to ours
				if(mComparisonLinks[j].appliesTo(port)){
					important=true;
					meas_pair=mComparisonLinks[j].measuredField();
					//meas_pair=mDataTable->otherColNameInComparisonLink(mComparisonLinks[j], outputname);
					break;
				}
			}
			if(important){
				if(srcname.size() && outputname.size() && meas_pair.size()){
					fprintf(ofile, "\t%s [shape=diamond];\n", outputname.c_str());
					fprintf(ofile, "\t%s [shape=diamond];\n", meas_pair.c_str());
					fprintf(ofile, "\t%s -> %s;\n", srcname.c_str(),outputname.c_str());
					fprintf(ofile, "\t%s:e -> %s:w [dir=both color=\"red:blue\"];\n", outputname.c_str(), meas_pair.c_str()); //style=dotted arrowhead=\"none\"
				}
			}
		}
	}
	
	fprintf(ofile,"}\n");
	fclose(ofile);
}

//---------------------------------------------------------------------------------------

#pragma mark Sensitivity analysis

void iWQModelLayout::localSensitivityAnalysis(double rel_deviance, std::string target, std::string filename)
{
	if(validity()<IWQ_VALID_FOR_RUN){
		printf("[Error]: Local sensitivity analysis failed: setup is not valid to run.\n");
		return;
	}
	
	//new columns
	std::vector<std::string> newcols;
	
	//basename
	std::string basename=target+std::string(" base");
	//newcols.push_back(basename);
	
	//make base run / with full warnings
	run();
	
	//copy results
	mDataTable->copyColumn(target, basename);
	
	//number of rows
	int ndata=mDataTable->numRows();
	
	//make a backup from parameters
	std::vector<double> par_backup=mCommonParameters->plainValues();
	std::vector<std::string> par_names=mCommonParameters->namesForPlainValues();
	
	//store outlets for scaling (later)
	double * baseval=mDataTable->portForColumn(basename);
	std::vector<double *> modvals;
	
	//make the sensitivity analysis
	std::vector<double> pars;
	int numpars=par_backup.size();
	for(int i=0; i<numpars; i++){
		//perturb a single parameter
		pars.assign(par_backup.begin(),par_backup.end());
		pars[i] *= 1.0 + rel_deviance;
		mCommonParameters->setPlainValues(pars);
		runmodel();
		std::stringstream s;
		s<<"SENSLOC_"<<target<<"_"<<par_names[i]<<"_"<<rel_deviance;
		std::string newname=s.str();
		mDataTable->copyColumn(target, newname);
		newcols.push_back(newname);
		modvals.push_back(mDataTable->portForColumn(newname));
	}
	
	//rescale the results to get relative sensitivity
	if(baseval){
		mDataTable->rewind();
		while(mDataTable->stepRow()!=-1){
			for(int i=0; i<numpars; i++){
				double * ptr=modvals[i];
				if(ptr){
					*ptr = (*ptr - *baseval) / *baseval / rel_deviance;
				}
			}
		}
	}
	else{
		printf("[Error]: Data error, only absolute sensitivity results were saved.\n");
	}
	
	// calculate correlation matrix between sensitivity functions
	//Eigen::MatrixXd corrmatrix (numpars, numpars);
	double ** corrmatrix=AllocMatrix(numpars);
		
	//get direct access to the new data columns
	std::vector<const std::vector<double>*> colvectors;
	const std::vector<double> * timevector;
	
	for(int i=0; i<numpars; i++){
		const std::vector<double> * ptr=mDataTable->vectorForColumn(newcols[i]);
		if(ptr){
			colvectors.push_back(ptr);
		}
		else{
			printf("[Error]: Data error, %s was skipped from the evaluation of sensitivity functions.\n",par_names[i].c_str());
		}
	}
	
	std::string tcol=mDataTable->timeColumn();
	timevector=mDataTable->vectorForColumn(tcol);
		
	for(int i=0; i<numpars; i++){
		for(int j=i; j<numpars; j++){
			double r=(i==j)?1.0:correlation(colvectors[i],colvectors[j]);
			corrmatrix[i][j]=r;
			if(i!=j){
				corrmatrix[j][i]=r;
			}
		}
	}
	
	std::vector<double> sensranks;
	for(int i=0; i<numpars; i++){
		sensranks.push_back(sqrt(sumsquares(colvectors[i])/(double)colvectors[i]->size()));
	}
	
	FILE * ofile=fopen(filename.c_str(),"w");
	if(!ofile){
		printf("[Error]: Failed to create %s.\n",filename.c_str());
	}
	else{
		fprintf(ofile,"LOCAL SENSITIVITY TEST for %s\nParamater perturbation=%d%%\n",target.c_str(),(int)(rel_deviance*100));
		fprintf(ofile,"\nSensitivity ranks:\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"%s\t%lf\n",par_names[i].c_str(),sensranks[i]);
		}
		fprintf(ofile,"\nCorrelation matrix between sensitivity functions:\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"\t%s",par_names[i].c_str());
		}
		fprintf(ofile,"\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"%s",par_names[i].c_str());
			for(int j=0; j<numpars; j++){
				fprintf(ofile,"\t%lf",corrmatrix[i][j]);
			}
			fprintf(ofile,"\n");
		}
		fprintf(ofile,"\nLocal sensitivity functions:\n");
		fprintf(ofile,"%s",tcol.c_str());
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"\t%s",par_names[i].c_str());
		}
		fprintf(ofile,"\n");
		for(int i=0; i<ndata; i++){
			if(timevector){
				fprintf(ofile,"%lf",timevector->at(i));
			}
			for(int j=0; j<numpars; j++){
				if(colvectors[j]){
					fprintf(ofile,"\t%lf",colvectors[j]->at(i));
				}
			}
			fprintf(ofile,"\n");
		}
		
		fclose(ofile);
	}
	
	
	//delete matrix
	DeleteMatrix(corrmatrix,numpars);
		
	//restore params & run results
	mCommonParameters->setPlainValues(par_backup);
	mDataTable->copyColumn(basename,target,false);
	
	//save results
	//mDataTable->writeToFile(filename,newcols);
	
	//cleanup new columns
	mDataTable->deleteColumn(basename);
	for(unsigned int i=0; i<newcols.size(); i++){
		mDataTable->deleteColumn(newcols[i]);
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::regionalSensitivityAnalysis(double rel_deviance, std::string target, std::string filename, int numsimulations)
{
	//storage for sim results
	std::vector<iWQVector> simruns;
	int numdata=mDataTable->numRows();
	
	//make a backup for parameters
	iWQVector par_backup=mCommonParameters->plainValues();
	std::vector<std::string> par_names=mCommonParameters->namesForPlainValues();	
	
	//basename
	std::string basename=target+std::string(" base");		//space is surely unique, it is not allowed in file headers
		
	//make base run / with full warnings
	run();
	
	//copy results
	mDataTable->copyColumn(target, basename);
	
	//accomodate parameter storage
	int numpars=par_backup.size();
	iWQVector pars;
	pars.assign(numpars,0.0);
	
	//ptr for the result column
	const iWQVector * resultcol=mDataTable->vectorForColumn(target);
	
	//make random parameters
	std::vector<iWQVector> randpars;
	iWQRandomLogNormalGenerator generator (0.0, 1.0);
	for(int i=0; i<numpars; i++){
		std::vector<double> parsamp;
		generator.setMean(par_backup[i]);
		generator.setStdev(rel_deviance*par_backup[i]);
		for(int j=0; j<numsimulations; j++){
			parsamp.push_back(generator.generate());
		}
		randpars.push_back(parsamp);
	}
	
	//make random simulations
	printf("Making random simulations...");
	
	int numfaulty=0;
	for(int k=0; k<numsimulations-numfaulty; k++){
		//feed parameters
		for(int i=0; i<numpars; i++){
			pars[i]=randpars[i][k];
		}
				
		//run simulation
		mCommonParameters->setPlainValues(pars);
		if(runmodel()){
			//stable solution
			printf(" %d",k+1+numfaulty);
		
			//copy results
			simruns.push_back(iWQVector(resultcol->begin(), resultcol->end()));
		}
		else{
			//unstable solution
			for(int i=0; i<numpars; i++){
				randpars[i].erase(randpars[i].begin()+k);	//get the faulty set out of the storage
			}
			k--;
			numfaulty++;
			printf(" ");
			int digits=(int)(log10(k+1+numfaulty));
			for(int l=0; l<=digits; l++){
				printf("-");
			}
		}
	}
	printf("\nReady\n");
	if(numfaulty){
		printf("%d numerically unstable solutions were omitted (~%d%%).\n",numfaulty,numfaulty*100/numsimulations);
	}
	
	numsimulations=simruns.size();	//actualize sim numbers
	
	//make parameter quantiles
	std::vector<iWQVector> quantiles;
	int nbins=(int)(sqrt(numsimulations))+1;
	for(int i=0; i<numpars; i++){
		std::vector<double> q;
		for(int j=0; j<nbins; j++){
			q.push_back(quantile(randpars[i], (j+0.5)/(double)nbins));
		}
		quantiles.push_back(q);
	}
	
	//analyse for each timestep
	iWQVector var_k;
	iWQVector curr_samp;
	std::vector<iWQVector> var_q;
	var_k.assign(numdata,0.0);
	var_q.assign(numpars,iWQVector());
	for(int i=0; i<numpars; i++){
		var_q[i].assign(numdata,0.0);
	}
	
	printf("Analysing variance...");
	curr_samp.assign(numsimulations,0.0);
	for(int i=0; i<numdata; i++){
		//prepare curr_samp
		for(int j=0; j<numsimulations; j++){
			curr_samp[j]=simruns[j][i];
		}
		
		//calculate variance
		var_k[i]=variance(&curr_samp);
				
		//calculate par dependent variance
		for(int j=0; j<numpars; j++){
			iWQVector fitted=loess(&randpars[j],&curr_samp,&quantiles[j]); //50/(double)numsimulations
			var_q[j][i]=variance(&fitted);			
		}		
	}
	printf("\tReady.\n");
	
	// calculate correlation matrix between sensitivity functions
	double ** corrmatrix=AllocMatrix(numpars);
		
	for(int i=0; i<numpars; i++){
		for(int j=i; j<numpars; j++){
			double r=(i==j)?1.0:correlation(&var_q[i],&var_q[j]);
			corrmatrix[i][j]=r;
			if(i!=j){
				corrmatrix[j][i]=r;
			}
		}
	}
	
	std::vector<double> sensranks;
	for(int i=0; i<numpars; i++){
		sensranks.push_back(sqrt(sumsquares(&var_q[i])/(double)var_q[i].size()));
	}
	
	//write out results
	FILE * ofile=fopen(filename.c_str(),"w");
	if(!ofile){
		printf("[Error]: Failed to create %s.\n",filename.c_str());
	}
	else{
		
		std::string tcol=mDataTable->timeColumn();
		const iWQVector * timevector=mDataTable->vectorForColumn(tcol);
		
		fprintf(ofile,"VARIANCE-BASED REGIONAL SENSITIVITY TEST for %s\nParameter sampling distribution: lognormal with %d%% stdev\n",target.c_str(),(int)(rel_deviance*100));
		
		fprintf(ofile,"\nSensitivity ranks:\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"%s\t%lf\n",par_names[i].c_str(),sensranks[i]);
		}
		fprintf(ofile,"\nCorrelation matrix between sensitivity functions:\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"\t%s",par_names[i].c_str());
		}
		fprintf(ofile,"\n");
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"%s",par_names[i].c_str());
			for(int j=0; j<numpars; j++){
				fprintf(ofile,"\t%lf",corrmatrix[i][j]);
			}
			fprintf(ofile,"\n");
		}
		
		fprintf(ofile,"\nRegional sensitivity functions:\n");
		fprintf(ofile,"%s\tVAR",tcol.c_str());
		for(int i=0; i<numpars; i++){
			fprintf(ofile,"\t%s",par_names[i].c_str());
		}
		fprintf(ofile,"\n");
		for(int i=0; i<numdata; i++){
			if(timevector){
				fprintf(ofile,"%lf",timevector->at(i));
			}
			fprintf(ofile,"\t%lf",var_k[i]);
			for(int j=0; j<numpars; j++){
				double val=(var_k[i]!=0.0)?sqrt(var_q[j][i]/var_k[i]):0.0;
				fprintf(ofile,"\t%lf",val);
			}
			fprintf(ofile,"\n");
		}
		fclose(ofile);
	}
	
	//delete matrix
	DeleteMatrix(corrmatrix,numpars);
	
	//restore params & run results
	mCommonParameters->setPlainValues(par_backup);
	mDataTable->copyColumn(basename,target,false); 
	mDataTable->deleteColumn(basename);
}

//---------------------------------------------------------------------------------------

#pragma mark Markov Chain Monte Carlo sampling of the evaluator function

//Utility to save the best parameters and the corresponding sample series

void iWQModelLayout::saveBestSolutionSoFar()
{
	std::string partempfile="_par_mcmc_best.txt";
	mCommonParameters->saveToFile(partempfile.c_str());
	
	//return;	//for debug purposes
	
	std::vector<std::string> seriessamples;
	for(int l=0; l<mEvaluatorMethods.size(); l++){
		std::vector<std::string> samples=mEvaluatorMethods[l]->sampleSeriesNames();
		seriessamples.insert(seriessamples.end(),samples.begin(), samples.end());
	}
	
	//get the series for the best solution so far
	std::map<std::string, std::vector<double> > sersampstorage;
	for(int o=0; o<mEvaluatorMethods.size(); o++){
		mEvaluatorMethods[o]->createSampleSeries(&sersampstorage);
	}
	//now save them as plain text files (numbers in a single column)
	for(int i=0; i<seriessamples.size(); i++){
		std::vector<double> data=sersampstorage[seriessamples[i]];
		std::string filename=seriessamples[i] + "_best.txt";
		FILE * fd=fopen(filename.c_str(),"w");
		int datasize=data.size();
		if(fd && datasize){
			for(int j=0; j<datasize; j++){
				fprintf(fd,"%lf\n",data[j]);
			}
			fclose(fd);
		}
	}
}

//---------------------------------------------------------------------------------------

std::string replaceParentheses(std::string s)
{
	size_t pos=s.find('[');
	if(pos!=std::string::npos){
		s[pos]='_';
	}
	pos=s.find(']');
	if(pos!=std::string::npos){
		s[pos]='_';
	}
	return s;
}

//Main MCMC routine

void iWQModelLayout::MCMC(int numrounds, int burnin, std::string filename, bool loadpropmatrix)
{
	if(!verify()){
		printf("[Error]: Model layout contains defects.\n");
		return;
	}
	
	printf("Markov-chain Monte Carlo experiment.\n");
		
	srand(time(0));
	int thinning=5;
	int nrounds=numrounds*thinning;
	int burn_in=burnin*thinning;
	
	int n=mCommonParameters->numberOfParams();
	double * parvals=new double [n];  //actual values
	double * nparvals=new double [n]; //storage for candidate values
	double * stdevs=new double [n];	  //sampler stdev
	std::vector<double> pars=mCommonParameters->plainValues();
	std::vector<double> orig_pars=pars;	//as a backup
	
	//to ease burning in
	printf("Initial optimization to ease burn-in.\n");
		
	mEvaluator->calibrate();
	
	//save best parameter values
	double besteval = mEvaluator->evaluate();
	saveBestSolutionSoFar();
			
	printf("Initial optimization finished.\nDoing MCMC:\n");	
	
	//start actual MCMC here
	double spreadfactor=1.0/12.0;
		
	printf("Thinning factor: %d\n",thinning);
	
	pars=mCommonParameters->plainValues();
	for(int i=0; i<n; i++){
		parvals[i]=pars[i];
		stdevs[i]=fabs(pars[i])*spreadfactor;	//was 1/12, arbitrary scaling for sample kernel TODO: replace
	}
	
	mEvaluator->printWarnings=false;
	
	double Pxt=mEvaluator->evaluate();
	double Pxi;
	FILE * ofile=fopen(filename.c_str(),"w");
	if(!ofile){
		printf("[Error]: Failed to open output file.\n");
		return;
	}
	std::vector<std::string> parnames=mCommonParameters->namesForPlainValues();
	
	fprintf(ofile,"step");
	for(int i=0; i<n; i++){
		std::string actparname=replaceParentheses(parnames[i]);	//remove [] from parnames because R does not like it
		fprintf(ofile,"\t%s",actparname.c_str());
	}
	fprintf(ofile,"\tEvaluation\talpha\tp\tburn_in\tacception\n");
	
	int thinindex=0;
	int accepted=0;
	int proposed=0;
	int adjustcycle=0;
	int faulty=0;
	double acception=0;
	
	//cached result output to HD
	std::string writecache="";
	char buf [1024];
	
	//BEGIN NEW
	int nsubsample = 1000;
	//storage for the subsample
	std::vector< std::vector<double> > subsample (n);
	for(int i=0; i<n; i++){
		subsample[i].assign(nsubsample, 0.0);
	}
	Eigen::MatrixXd SIGMA;
	Eigen::MatrixXd L_SIGMA;
	Eigen::MatrixXd SIGMA_OLD;
	Eigen::MatrixXd L_SIGMA_OLD;
	Eigen::MatrixXd SIGMA_ALT;
	
	bool initblankrun = true;
	double blankrunlimit = 0.15;
	
	std::string propmatname = "_proposal_matrix.txt";
	
	if(loadpropmatrix){
		SIGMA=LoadMatrix(propmatname);
		if(SIGMA.cols()>0 && SIGMA.rows()>0){
			try{
				L_SIGMA = choleskyDecomposition(SIGMA);
				initblankrun=false;
			}
			catch(...){
				printError("Cholesky decomposition failed, not using imported proposal matrix.",NULL);
				initblankrun=true;;
			}
			if(!initblankrun){
				std::cout<<"*** Using imported proposal matrix ***\n"<<std::endl;
			}
		}
	}
	
	std::vector<double> mus (n);	//expected values
	std::vector<double> zs (n);		//i.i.d. std normals
		
	iWQRandomNormalGenerator N;
	
	double c = 1.0;
	double maxr2 = 0.95;
	
	//END NEW
	
	for(int i=0; i<nrounds; i++){
		//generate a new candidate
		//samplerKernel(n, parvals, nparvals, stdevs);
		
		//BEGIN NEW
		if(initblankrun && acception<blankrunlimit){
			//blank draw in 1st period
			
			samplerKernel(n, parvals, nparvals, stdevs);
			
		}
		else{
			//proper multivariate draw
			for(int j=0; j<n; j++){
				mus[j]=parvals[j];
				zs[j]=N.generate();
			}
			std::vector<double> newpars = multivariateNormal(L_SIGMA, zs, mus);
			for(int j=0; j<n; j++){
				nparvals[j]=newpars[j];
			}
		}
		//END NEW
		
		//evaluate it
		mCommonParameters->setPlainValues(nparvals,n);
		Pxi=mEvaluator->evaluate();
		
		if(Pxi!=DBL_MAX){ //if valid run
			//first let's see its absolute performace
			if(Pxi<besteval){
				saveBestSolutionSoFar();
				besteval=Pxi;
			}
			//decide on acceptance
			proposed++;
			double a;
			if(!mEvaluatorMethods[0]->isLogScale()){
				a=(Pxi!=0.0 && !(std::isnan(Pxi) || std::isinf(Pxi)))?(Pxt/Pxi):-1.0; 	//reversed, since we want to minimize	//was pow(Pxt/Pxi,n), but why?
				if(std::isnan(a) || std::isinf(a)){
					printf("\nError in linear draw: a=%lf Pxt=%lf Pxi=%lf n=%d\n",a,Pxt,Pxi,n);
				}
			}
			else{
				a=(Pxi>Pxt)?exp((Pxt - Pxi)):1.0; 	//reversed, since we want to minimize	//was (Pxt-Pxi)*((double)n), but why?
				if(std::isnan(a) || std::isinf(a)){
					printf("\nError in log draw: a=%lf Pxt=%lf Pxi=%lf n=%d\n",a,Pxt,Pxi,n);
				}
			}
			double p=urand();
			if(p<a){				
				//copy if needed
				Pxt=Pxi;
				for(int j=0; j<n; j++){
					parvals[j]=nparvals[j];
				}
				accepted++;
			}
			
			//BEGIN NEW
			//save the final values to the sample container
			if(initblankrun || i<burn_in){
				for(int j=0; j<n; j++){
					subsample[j][proposed-1]=parvals[j];
				}
			}
			//END NEW
			
			thinindex++;
			if(thinindex>=thinning){
				//save result
				memset(buf, '\0', sizeof(buf) );
				sprintf(buf,"%d",i/thinning);
				writecache+=buf;
				for(int j=0; j<n; j++){
					memset(buf, '\0', sizeof(buf) );
					sprintf(buf,"\t%g",parvals[j]);
					writecache+=buf;
				}
				memset(buf, '\0', sizeof(buf) );
				sprintf(buf,"\t%lf\t%lf\t%lf\t%c\t%d%%\n",Pxt,a,p,(i>=burn_in?'1':'0'),(int)(acception*100));
				writecache+=buf;
				
				if(initblankrun){
					printf("*");
				}
				if(i<burn_in){
					printf("*");
				}
				printf("%d ",i/thinning);
				fflush(stdout);
				thinindex=0;				
			}
		}
		else{
			//proposed++;		//need to consider numerically unstable proposals too
			faulty++;
			int digits=1;
			if(i>0){ 
				digits=(int)log10(i)+1;
			}
			for(int j=0; j<digits; j++){
				printf("-");
			}
			printf(" ");
			fflush(stdout);
			i--;
		}
		
		//rotate acception stats
		if(proposed>=nsubsample){
			acception=accepted/(double)proposed;

			//review the initial period
			if(initblankrun && acception>blankrunlimit){
				initblankrun=false;
				//make the burn_in phase longer : not good for post processing
				//announce this
				printf("\n\n*** INITIAL SCALING PERIOD IS OVER. ***\n\n");
				
				//here fix the covar matrix and its initial scale
				try{
					SIGMA = covarMatrix2(subsample, maxr2);
					if(!isfinite(SIGMA)){
						//invalid covariance matrix as start
						//initblankrun=true;
						//proposed=0;   //?
						printError("Invalid covariance matrix, restarting initial sampling.",NULL);
						std::cout << SIGMA <<std::endl;
					}
				}
				catch(...){
					printf("\nError happened during the construction of the initial covariance matrix. Process aborted.\n");
					initblankrun=true;
					return;
				}
				c = 1.0; 	//restart the scaling
			}
			
			if(i<burn_in){
				if(initblankrun){
					//initial scaling phase
					double std_mod=1.0;
					if(acception<0.15){		// Limits are taken from Gelman, Roberts and Gilks (1996)
						std_mod=0.8;		// Before this we used 0.3 and 0.7 to set the goal to 50% 
					}
					else if(acception>0.4){
						std_mod=1.2;
					}
					c *= std_mod;
					printf("\n*** Adjusting initial proposal width (to %lf) ***\n",c);
					
					//recalculate stdevs
					for(int j=0; j<n; j++){
						stdevs[j]*=std_mod;
					}
				}
				else{
					//if not initblankrun
					//more sophisticated proposal distribution tuning
					//covar matrix update
					double updateratio = 0.1;
					
					SIGMA_OLD = SIGMA;
					try{
						Eigen::MatrixXd SIGMA_NEW = covarMatrix2(subsample, maxr2);
						//filter if SIGMA_NEW contains anything illegal
						if(isfinite(SIGMA_NEW)){
							SIGMA = ((1.0 - updateratio) * SIGMA) + (updateratio * SIGMA_NEW);	//"dampened sequential update", does not do anything in a converged state
						}else{
							printError("Covariance matrix calculation failed (1), reverting to backup matrix.",NULL);
						}
					}
					catch(...){
						printError("Covariance matrix calculation failed (2), reverting to backup matrix.",NULL);
						SIGMA = SIGMA_OLD;
					}
					
					if(acception<0.15){		
						c *= 0.9;
						printf("\n*** Narrowing proposal distribution (to %lf) ***\n",c);
					}
					else if(acception>0.4){
						c *= 1.1;
						printf("\n*** Widening proposal distribution (to %lf) ***\n",c);
					}
					
					SIGMA_ALT = (c*c) * SIGMA;
					L_SIGMA_OLD = L_SIGMA;
					try{
						L_SIGMA = choleskyDecomposition(SIGMA_ALT);
					}
					catch(...){
						printError("Cholesky decomposition failed (1), reverting to backup matrix.",NULL);
						L_SIGMA = L_SIGMA_OLD;
					}
					if(!isfinite(L_SIGMA)){
						printError("Cholesky decomposition failed (2), reverting to backup matrix.",NULL);
						L_SIGMA = L_SIGMA_OLD;
					}
					//empty subsamples
					//not necessary
					//DEBUG: occasional instability is caused by extreme low values for certain parameters
					//SaveMatrix(SIGMA_OLD, "_original_proposal_matrix.txt");
					//SaveMatrix(SIGMA, "_updated_proposal_matrix.txt");
					//END DEBUG
				}
			}
			
			printf("\n*** Last %d rounds acception statistics: %d%% ***\n",nsubsample,(int)(acception*100));
			printf("*** Current best likelihood point: [%lf] ***\n",Pxt);
			//restart evaluation
			proposed=0;
			accepted=0;
			faulty=0;
		}
		
		if(i<burn_in && (proposed==0 && faulty>=nsubsample)){
			//deadlock: no numerically valid proposals within a whole period
			printf("\n*** DEADLOCK REACHED WITH ADAPTING THE PROPOSAL DISTRIBUTION ***\n");
			//printf("Reverting to the last known uncorrelated proposal distribution!\n");
			printf("Actual parameter values:\n");
			std::vector<double> par_vals=mCommonParameters->plainValues();
			std::vector<std::string> par_names=mCommonParameters->namesForPlainValues();
			for(int i=0; i<par_names.size(); i++){
				printf("\t%s:\t%lf\n",par_names[i].c_str(),par_vals[i]);
			}
			//blankrunlimit=1.0; //never reached
			initblankrun=true;
			faulty=0;
			return;
		}
		
		//rotate write cache
		if(writecache.size()>=10240){	//write out in 10KB chunks
			fprintf(ofile,"%s",writecache.c_str());
			fflush(ofile);
			writecache="";
		}
	}//end of iteration
	
	//flush remaining write cache (if any)
	if(writecache.size()){
		fprintf(ofile,"%s",writecache.c_str());
	}
	
	//save proposal distribution to file
	SaveMatrix(SIGMA, propmatname);
	
	fclose(ofile);
	
	delete [] parvals;
	delete [] nparvals;
	delete [] stdevs;
	
	//restore original parameter values
	mCommonParameters->setPlainValues(orig_pars);
	mEvaluator->printWarnings=true;
	//mSolver->setMinStepLength(orig_minstep);
	
	//BEGIN NEW: close series samples files		//WAS REMOVED TO SPEED UP MCMC
	//close the series sample files
	/*std::map<std::string, FILE *>::iterator it;
	for(it=seriessamplefiles.begin(); it!=seriessamplefiles.end(); ++it){
		int endtag=0;
		fwrite(&endtag,sizeof(int),1,it->second);
		fclose(it->second);
	}*/
	//END NEW: close series samples files
	
	printf("\nMCMC sampling finished.\n");
	
	return;
}

//Main MCMC routine

void iWQModelLayout::MCMC_Haario(int numrounds, int burnin, std::string filename)
{
	if(!verify()){
		printf("[Error]: Model layout contains defects.\n");
		return;
	}
	
	printf("Markov-chain Monte Carlo experiment (Haario\'s algorithm).\n");
		
	srand(time(0));
	int thinning=5;
	int nrounds=numrounds*thinning;
	int burn_in=burnin*thinning;
	
	int n=mCommonParameters->numberOfParams();
	double * parvals=new double [n];  //actual values
	double * nparvals=new double [n]; //storage for candidate values
	double * stdevs=new double [n];	  //sampler stdev
	std::vector<double> pars=mCommonParameters->plainValues();
	std::vector<double> orig_pars=pars;	//as a backup
	
	//to ease burning in
	printf("Initial optimization to ease burn-in.\n");
	mEvaluator->calibrate();
	
	//save best parameter values
	double besteval = mEvaluator->evaluate();
	saveBestSolutionSoFar();
	
	printf("Initial optimization finished.\nDoing MCMC (Haario):\n");	
	
	//start actual MCMC here
	double spreadfactor=1.0/12.0;
		
	printf("Thinning factor: %d\n",thinning);
	
	pars=mCommonParameters->plainValues();
	for(int i=0; i<n; i++){
		parvals[i]=pars[i];
		stdevs[i]=fabs(pars[i])*spreadfactor;
	}
	
	mEvaluator->printWarnings=false;
		
	double Pxt=mEvaluator->evaluate();
	double Pxi;
	FILE * ofile=fopen(filename.c_str(),"w");
	if(!ofile){
		printf("[Error]: Failed to open output file.\n");
		return;
	}
	std::vector<std::string> parnames=mCommonParameters->namesForPlainValues();
	
	fprintf(ofile,"step");
	for(int i=0; i<n; i++){
		std::string actparname=replaceParentheses(parnames[i]);	//remove [] from parnames because R does not like it
		fprintf(ofile,"\t%s",actparname.c_str());
	}
	fprintf(ofile,"\tEvaluation\talpha\tp\tburn_in\tacception\n");
	
	int thinindex=0;
	int accepted=0;
	int proposed=0;
	int adjustcycle=0;
	double acception=0;
	
	//cached result output to HD
	std::string writecache="";
	char buf [1024];
	
	//storage for the parameter trace
	std::vector< std::vector<double> > subsample (n); //subsample[p][j] where p is parameter, j is iteration index
	
	Eigen::MatrixXd SIGMA;
	Eigen::MatrixXd L_SIGMA;
	Eigen::MatrixXd SIGMA_ALT;
	
	bool initblankrun = true;
	double blankrunlimit = 0.25;
	
	std::vector<double> mus (n);	//expected values
	std::vector<double> zs (n);		//i.i.d. std normals
		
	iWQRandomNormalGenerator N;
	
	double c = 1.0;
	
	for(int i=0; i<nrounds; i++){
		//generate a new candidate
		if(initblankrun && acception<blankrunlimit){
			//blank draw in 1st period
			samplerKernel(n, parvals, nparvals, stdevs);
		}
		else{
			//proper multivariate draw
			for(int j=0; j<n; j++){
				mus[j]=parvals[j];
				zs[j]=N.generate();
			}
			std::vector<double> newpars = multivariateNormal(L_SIGMA, zs, mus);
			for(int j=0; j<n; j++){
				nparvals[j]=newpars[j];
			}
		}
				
		//evaluate it
		mCommonParameters->setPlainValues(nparvals,n);
		Pxi=mEvaluator->evaluate();
		
		if(Pxi!=DBL_MAX){ //if valid run
			//first let's see its absolute performace
			if(Pxi<besteval){
				saveBestSolutionSoFar();
				besteval=Pxi;
			}
			//decide on acceptance
			proposed++;
			double a;
			if(!mEvaluatorMethods[0]->isLogScale()){
				a=(Pxi!=0.0 && !(std::isnan(Pxi) || std::isinf(Pxi)))?(Pxt/Pxi):-1.0; 	//reversed, since we want to minimize	//was pow(Pxt/Pxi,n), but why?
				if(std::isnan(a) || std::isinf(a)){
					printf("\nError in linear draw: a=%lf Pxt=%lf Pxi=%lf n=%d\n",a,Pxt,Pxi,n);
				}
			}
			else{
				a=(Pxi>Pxt)?exp((Pxt - Pxi)):1.0; 	//reversed, since we want to minimize	//was (Pxt-Pxi)*((double)n), but why?
				if(std::isnan(a) || std::isinf(a)){
					printf("\nError in log draw: a=%lf Pxt=%lf Pxi=%lf n=%d\n",a,Pxt,Pxi,n);
				}
			}
			double p=urand();
			if(p<a){				
				//copy if needed
				Pxt=Pxi;
				for(int j=0; j<n; j++){
					parvals[j]=nparvals[j];
				}
				accepted++;
			}
			
			//save the final values to the sample container
			if(p<a && i<0.75 * nrounds){
				for(int j=0; j<n; j++){
					subsample[j].push_back(parvals[j]);
				}
			}
			
			thinindex++;
			if(thinindex>=thinning){
				//save result
				memset(buf, '\0', sizeof(buf) );
				sprintf(buf,"%d",i/thinning);
				writecache+=buf;
				for(int j=0; j<n; j++){
					memset(buf, '\0', sizeof(buf) );
					sprintf(buf,"\t%g",parvals[j]);
					writecache+=buf;
				}
				memset(buf, '\0', sizeof(buf) );
				sprintf(buf,"\t%lf\t%lf\t%lf\t%c\t%d%%\n",Pxt,a,p,(i>=burn_in?'1':'0'),(int)(acception*100));
				writecache+=buf;
				
				if(initblankrun){
					printf("*");
				}
				if(i<burn_in){
					printf("*");
				}
				printf("%d ",i/thinning);
				fflush(stdout);
				thinindex=0;				
			}
		}
		else{
			//proposed++;		//need to consider numerically unstable proposals too
			int digits=1;
			if(i>0){ 
				digits=(int)log10(i);
			}
			for(int j=0; j<digits; j++){
				printf("-");
			}
			printf(" ");
			i--;
		}
		
		//rotate acception stats
		int rotatefreq = 100;
		if(proposed>=rotatefreq){				//update after 100 proposals
			acception=accepted/(double)proposed;

			//review the initial period
			if(initblankrun && acception>blankrunlimit){
				initblankrun=false;
				//announce this
				printf("\n\n*** INITIAL SCALING PERIOD IS OVER. ***\n\n");
				
				//here setup the initial covar matrix
				int ninitsample = (n>0)?subsample[0].size():1;
				int ninitcovar = (ninitsample>200)?200:ninitsample;	     //use only the last 200 items from the subsample 
				printf("Calculating the covariance matrix from the last %d accepted sets.\n",ninitsample);
				SIGMA = covarMatrix(subsample, ninitcovar);
				printf("Ready.\n");
				//discard the existing subsample (to exclude the uncorrelated part from future samples)
				/*for(int j=0; j<n; j++){
					subsample[j].clear();
				}*/
				c = 1.0; 	//restart the scaling too
			}
			
			if(i<burn_in){
				if(initblankrun){
					//initial scaling phase
					double std_mod=1.0;
					if(acception<0.2){		// Limits are taken from Gelman, Roberts and Gilks (1996)
						std_mod=0.9;		// Before this we used 0.3 and 0.7 to set the goal to 50% 
					}
					else if(acception>0.4){
						std_mod=1.1;
					}
					c *= std_mod;
					printf("\n*** Adjusting initial proposal width (to %lf) ***\n",c);
					
					//recalculate stdevs
					for(int j=0; j<n; j++){
						stdevs[j]*=std_mod;
					}
				}
				else{
					//if not initblankrun
					//more sophisticated proposal distribution tuning
					//covar matrix update
					int N = (n>0)?subsample[0].size():0;
					if(N>0){
						double correctratio = (N - 1.0) / (double)N;
						SIGMA = correctratio * covarMatrix(subsample, (int)(N*0.5));	//continuous updating of the covariance matrix
					}					
					c=1.0;
					if(acception<0.15){		
						c = 0.8;
						printf("\n*** Temporarily narrowing proposal distribution scale (to %lf) ***\n",c);
					}
					else if(acception>0.4){
						c = 1.2;
						printf("\n*** Temporarily widening proposal distribution scale (to %lf) ***\n",c);
					}
					
					SIGMA_ALT = c * SIGMA;
					
					L_SIGMA = choleskyDecomposition(SIGMA_ALT);
				}
			}
			
			printf("\n*** Last %d rounds acception statistics: %d%% ***\n",rotatefreq,(int)(acception*100));
			printf("*** Current best likelihood point: [%lf] ***\n",Pxt);
			//restart evaluation
			proposed=0;
			accepted=0;
		}
		
		//rotate write cache
		if(writecache.size()>=10240){	//write out in 10KB chunks
			fprintf(ofile,"%s",writecache.c_str());
			fflush(ofile);
			writecache="";
		}
	}
	
	//flush remaining write cache (if any)
	if(writecache.size()){
		fprintf(ofile,"%s",writecache.c_str());
	}
	
	fclose(ofile);
	
	delete [] parvals;
	delete [] nparvals;
	delete [] stdevs;
	
	//restore original parameter values
	mCommonParameters->setPlainValues(orig_pars);
	mEvaluator->printWarnings=true;
		
	printf("\nMCMC_HAARIO sampling finished.\n");
	
	return;
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::runOnSample(std::string samplefilename, std::string outputfilename)
{
	//open the sample file
	iWQDataTable * datatable=new iWQDataTable (samplefilename);
		
	//check if we have anything in the sample
	if(datatable->numRows()<=0){
		printf("[Error]: Failed to load parameter sample from file \"%s\".\n",samplefilename.c_str());
		delete datatable;
		return;
	}
	
	datatable->addColumn("Numerical_stability");
	double * qualityflagptr=datatable->portForColumn("Numerical_stability");
	
	//wire in parameters
	std::vector<std::string> paramnames=mCommonParameters->namesForPlainValues();
		
	int nrows=datatable->numRows();
	int nparams=paramnames.size();
	
	double ** paramloc=new double * [nparams];
	double * paramvalues=new double [nparams];
	
	for(int i=0; i<nparams; i++){
		double * loc=datatable->portForColumn(replaceParentheses(paramnames[i]));
		double * alt_loc = NULL;
		if(paramnames[i].compare(replaceParentheses(paramnames[i]))!=0){
			alt_loc=datatable->portForColumn(paramnames[i]);
		}
		if(!loc && !alt_loc){
			printf("[Error]: Parameter \"%s\" was not found in the sample file.\n",paramnames[i].c_str());
			delete [] paramloc;
			delete [] paramvalues;
			delete datatable;
			return;
		}
		if(loc && alt_loc){
			printf("[Error]: Parameter \"%s\" was found both with [XX] and _XX_ syntax.\n",paramnames[i].c_str());
			delete [] paramloc;
			delete [] paramvalues;
			delete datatable;
			return;
		}
		if(loc){
			paramloc[i]=loc;
		}
		else{
			paramloc[i]=alt_loc;
		}
	}
	
	int numfaulty=0;
	printf("Running sample simulations...");
	for(int r=0; r<nrows; r++){
		datatable->setRow(r);
		//set the new parameters
		for(int i=0; i<nparams; i++){
				paramvalues[i] = *paramloc[i];
		}
		mCommonParameters->setPlainValues(paramvalues, nparams);
		
		if(mSeriesInterface){
			mSeriesInterface->refreshInputs();
		}
		
		//run and record the results 
		if(runmodel()){
			//stable solution
			printf(" %d",r);
			*qualityflagptr=1;
		}
		else{
			//unstable solution: should not happen
			numfaulty++;
			printf(" ");
			int digits=(int)log10(r+1);
			for(int l=0; l<=digits; l++){
				printf("-");
			}
			*qualityflagptr=0;
		}
		
		if(mSeriesInterface){
			mSeriesInterface->refreshOutputs();
		}
		//step to the next parameter sample row
		datatable->stepRow();
	}
	
	datatable->writeToFile(outputfilename);
	
	printf("\nReady\n");
	
	delete datatable;
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::createBestSeries(std::string parbestfilename)
{
	//if we have a best parameter file
	if(parbestfilename.size()){
		mCommonParameters->initFromFile(parbestfilename);
		//run the model
		double eval = mEvaluator->evaluate();
		saveBestSolutionSoFar();		
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::runStandardSeriesOnSample(std::string samplefilename, int desiredrowcount, bool predictivemode, bool binary)
{	
	//open the sample file
	iWQDataTable * datatable=new iWQDataTable (samplefilename);
		
	//check if we have anything in the sample
	if(datatable->numRows()<=0){
		printf("[Error]: Failed to load parameter sample from file \"%s\".\n",samplefilename.c_str());
		delete datatable;
		return;
	}
	
	double * burninptr = datatable->portForColumn("burn_in");
	
	//furnish series sample files
	std::vector<std::string> seriessamples;
	for(int l=0; l<mEvaluatorMethods.size(); l++){
		std::vector<std::string> samples=mEvaluatorMethods[l]->sampleSeriesNames();
		seriessamples.insert(seriessamples.end(),samples.begin(), samples.end());
	}
	std::map<std::string, FILE *> seriessamplefiles;
	for(int k=0; k<seriessamples.size(); k++){
		std::string fname;
		FILE * fd;
		if(binary){
			fname="series_"+seriessamples[k] + ".series";
			fd=fopen(fname.c_str(),"wb");
		}
		else{
			fname="series_"+seriessamples[k] + ".txt";
			fd=fopen(fname.c_str(),"w");
		}
		if(!fd){
			printf("[Error]: Failed to open %s.\n",fname.c_str());
		}
		else{
			seriessamplefiles[seriessamples[k]]=fd;
		}
	}
	
	//get the number of burned-in rows
	int thinning=1;
	if(burninptr && desiredrowcount>0){
		int numrows = datatable->numRows();
		int binrows = 0;
		for(int i=0; i<numrows; i++){
			datatable->setRow(i);
			if(*burninptr>0){
				binrows++;
			}
		}
		thinning = binrows / desiredrowcount;
		if(thinning<1){
			thinning=1;
		}
	}
	
	//wire in parameters
	std::vector<std::string> paramnames=mCommonParameters->namesForPlainValues();
		
	int nrows=datatable->numRows();
	int nparams=paramnames.size();
	
	double ** paramloc=new double * [nparams];
	double * paramvalues=new double [nparams];
	
	for(int i=0; i<nparams; i++){
		std::string rstyle_parname = replaceParentheses(paramnames[i]);
		double * loc=datatable->portForColumn(rstyle_parname);
		double * alt_loc = NULL;
		if(paramnames[i].compare(rstyle_parname)!=0){
			alt_loc=datatable->portForColumn(paramnames[i]);
		}
		if(!loc & !alt_loc){
			printf("[Error]: Parameter \"%s\"/\"%s\" was not found in the sample file.\n",rstyle_parname.c_str(),paramnames[i].c_str());
			delete [] paramloc;
			delete [] paramvalues;
			delete datatable;
			return;
		}
		if(loc && alt_loc){
			printf("[Error]: Parameter \"%s\" was also found as \"%s\".\n",paramnames[i].c_str(), rstyle_parname.c_str());
			delete [] paramloc;
			delete [] paramvalues;
			delete datatable;
			return;
		}
		if(loc){
			paramloc[i]=loc;
		}
		else{
			paramloc[i]=alt_loc;
		}
	}
	
	//set predictive mode
	for(int i=0; i<mComparisonLinks.size(); i++){
		mComparisonLinks[i].setPredictiveMode(predictivemode);
	}
	mEvaluator->setPredictiveMode(predictivemode);
	
	int numfaulty=0;
	printf("Running sample simulations (purely predictive mode)...\n");
	mEvaluator->printWarnings=false;
	for(int r=0; r<nrows; r++){
		datatable->setRow(r);
		//printf("Scanning row #%d\n",r);
		if(burninptr && *burninptr>0){
			//set the new parameters
			for(int i=0; i<nparams; i++){
					paramvalues[i] = *paramloc[i];
			}
			mCommonParameters->setPlainValues(paramvalues, nparams);
			
			if(mSeriesInterface){
				mSeriesInterface->refreshInputs();
			}
			
			//run and record the results 
			//printf("\tEvaluating row #%d\n",r);
			double evalres = mEvaluator->evaluate();
			//printf("Row #%d evaluated to %lf\n",r,evalres);
			if(evalres!=DBL_MAX){
			//if(runmodel()){
				//stable solution
				printf(" %d",r);
				/*if(evalres==DBL_MAX){
					printf("!");
				}*/
				fflush(stdout);
				
				//since we have an output, call the series sampler and write out its results
				std::map<std::string, std::vector<double> > sersampstorage;
				for(int o=0; o<mEvaluatorMethods.size(); o++){
					mEvaluatorMethods[o]->createSampleSeries(&sersampstorage);
				}
				for(int o=0; o<seriessamples.size(); o++){
					std::map<std::string, std::vector<double> >::iterator itsamp=sersampstorage.find(seriessamples[o]);
					//std::vector<double> data=sersampstorage[seriessamples[o]];	//prevent copying the whole vector
					FILE * fd=seriessamplefiles[seriessamples[o]];
					if(itsamp!=sersampstorage.end()){
						if(fd && itsamp->second.size()){
							//binary output
							int datasize=itsamp->second.size(); //data.size();
							if(binary){
								fwrite(&datasize,1,sizeof(int),fd);
								for(int m=0; m<datasize; m++){
									float val=itsamp->second[m]; //data[m];		//single precision output
									fwrite(&val,1,sizeof(float),fd);
								}
							}
							else{
								for(int m=0; m<datasize; m++){
									float val=itsamp->second[m]; //data[m];		//single precision output
									if(m>0){
										fprintf(fd,"\t");
									}
									fprintf(fd, "%f", val);
								}
								fprintf(fd,"\n");
							}
						}
						else{
							printf("\n[Error]: Corrupt series sample: FD=%p, data size=%zd\n",fd, itsamp->second.size());
						}
					}
					else{
						printf("\n[Error]: Could not find series sample for %s.\n",seriessamples[o].c_str());
					}
				}
				
				if(mSeriesInterface){
					mSeriesInterface->refreshOutputs();
				}
			}
			else{
				//unstable solution: should not happen
				numfaulty++;
				printf(" ");
				int digits=(int)log10(r+1);
				for(int l=0; l<=digits; l++){
					printf("-");
				}
				//DEBUG: spit out parameters
			}
			r+=(thinning-1);
		}
	}
	mEvaluator->printWarnings=true;
	
	//close the series sample files
	std::map<std::string, FILE *>::iterator it;
	for(it=seriessamplefiles.begin(); it!=seriessamplefiles.end(); ++it){
		int endtag=0;
		fwrite(&endtag,sizeof(int),1,it->second);
		fclose(it->second);
	}
	
	//restore predictive mode to OFF (default)
	for(int i=0; i<mComparisonLinks.size(); i++){
		mComparisonLinks[i].setPredictiveMode(false);
	}
	mEvaluator->setPredictiveMode(false);
	
	printf("\nReady\n");
	
	delete datatable;
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::furnishUNCSIM(std::string dirname)
{
	if(validity()<IWQ_VALID_FOR_CALIBRATE){
		printf("[Error]: Cannot make UNCSIM configuration files, because the specified layout is not valid for calibration.\n");
		return;
	}
	
	if(mkdir(dirname.c_str(),777)==-1){
		//failed to create
		//check if it exists
		struct stat stFileinfo;
		if(stat(dirname.c_str(), &stFileinfo)!=0){
			printf("[Error]: No access to the output directory.\n");
			return;
		}
		else{
			//we have info, so it exists
		}
	}
	
	// Start making the files
	std::string configfilename=dirname+std::string("/config.txt");
	std::string pardeffilename=dirname+std::string("/pardef.txt");
	std::string likelideffilename=dirname+std::string("/likelidef.txt");
	
	FILE * ofile=NULL;
	
	//configuration file
	ofile=fopen(configfilename.c_str(), "w");
	if(ofile){
		fprintf(ofile, 	"Model\tExternal\n"
						"External_ModelInFile\tuncsim-model.in\n"
						"External_ModelOutFile\tuncsim-model.out\n"
						"External_ModelExecFile\tuncsim-model.bat\n"
						"\n"
						"MaxIter\t5000\n"
						"TransposeQuant\tT\n"
						"\n"
						"ParDefFile\tpardef.txt\n"
						"LikeliDefFile\tlikelidef.txt\n"
						"\n"
						"ResValFile\tout_resval.txt\n"
						"TSResValFile\tout_tsresval.txt\n"
						"ParDefOutFile\tout_pardef.txt\n"
						"ResidValFile\tout_residval.txt\n"
						"ParTraceFile\tout_partrace.txt\n");
		fclose(ofile);
		ofile=NULL;
	}
	else{
		printf("[Error]: Failed to create the configuration file.\n");
	}
	
	//parameter definition file
	ofile=fopen(pardeffilename.c_str(), "w");
	if(ofile){
		fprintf(ofile, 	"Name\tValue\tMinimum\tMaximum\tScale\tUncRange\tIncrement\tActSens\tActEstim\tUnit\tDescription\n");
		int numpars=mCommonParameters->numberOfParams();
		std::vector<double> values=mCommonParameters->plainValues();
		std::vector<std::string> names=mCommonParameters->namesForPlainValues();
		
		//actual parameters
		for(int i=0; i<numpars; i++){
			iWQLimits lim;
			lim.min=0.0;
			lim.max=values[i]*10.0;
			//printf("Querying %s...\n",names[i].c_str());
			if(mCommonParameters->hasLimitsForParam(names[i])){
				lim=mCommonParameters->limitsForParam(names[i]);
				//printf("%s has limits: %lf, %lf\n",names[i].c_str(),lim.min, lim.max);
			}
			fprintf(ofile,"%s\t%lf\t%lf\t%lf\t1\t%lf\t%lf\tT\tT\t-\t\n",names[i].c_str(), values[i], lim.min, lim.max, values[i], (lim.max-lim.min)/1000.0);
		}
		
		//likelihood parameters
		fprintf(ofile,	"_sd\t0.5\t0.5\t0.5\t0.5\t0\t1\tF\tF\t-\n"
						"_l1\t0.3\t0.3\t0.3\t1\t0\t1\tF\tF\t-\n"
						"_l2\t0.1\t0.1\t0.1\t1\t0\t1\tF\tF\t-\n");
		
		fclose(ofile);
		ofile=NULL;
	}
	else{
		printf("[Error]: Failed to create parameter definition file.\n");
	}
	
	//likelihood definition file
	ofile=fopen(likelideffilename.c_str(), "w");
	if(ofile){
		fprintf(ofile, "ResCode\tData\tTransformation\tTransPar1\tTransPar2\tDistribution\tDistPar1\tDistPar2\n");
		//print out measured data from the comparison links
		for(int i=0; i<mComparisonLinks.size(); i++){
			//get colnames
			std::vector<std::string> data=mDataTable->UNCSIMdata(mComparisonLinks[i].measuredField(), mComparisonLinks[i].modelField());
			for(int j=0; j<data.size(); j++){
				fprintf(ofile, "%s\tBoxCox\t_l1\t_l2\tNormal\t0\t_sd\n", data[j].c_str());
			}
		}
		fclose(ofile);
	}
}

//---------------------------------------------------------------------------------------

void iWQModelLayout::printModelInfo(std::string name)
{
	//decide if this is a specific unit ID
	iWQModel * m=NULL;
	for(int i=0; i<mModels.size(); i++){
		if(mModels[i]){
			std::string act_id=mModels[i]->modelId();
			if(act_id==name){
				m=mModels[i];
				break;
			}
		}
	}
	//if not, try to allocate it
	bool demo=false;	//will show if we need the model after the query
	if(!m){
		if(mModelFactory){
			m=mModelFactory->newModelOfType(name);
			if(m){
				demo=true;
			}
		}
	}
	
	//check if we have something to ask
	if(!m){
		printf("[Error]: \"%s\" is neither a unit ID nor a type name.\n",name.c_str());
		return;
	}
	
	//ask the model
	if(demo){
		printf("*** Properties of type \"%s\" ***\n",name.c_str());
	}
	else{
		printf("*** Properties of unit \"%s\" ***\n",name.c_str());
		printf("\tType:        %s\n",m->modelType().c_str());
	}
	iWQStrings inputs, outputs, vars, params;
	inputs=m->inputDataHeaders();
	outputs=m->outputDataHeaders();
	params=m->parameters();
	int numvariables=m->numVariables();
	
	printf("\tInputs:      ");
	for(int i=0; i<inputs.size(); i++){
		if(i>0){
			printf(", ");
		}
		printf("%s",inputs[i].c_str());
	}
	printf("\n");
	
	printf("\tVariables:   ");
	for(int i=0; i<numvariables; i++){
		if(i>0){
			printf(", ");
		}
		printf("%s",outputs[i].c_str());
	}
	printf("\n");
	
	printf("\tFluxes:      ");
	for(int i=numvariables; i<outputs.size(); i++){
		if(i>numvariables){
			printf(", ");
		}
		printf("%s",outputs[i].c_str());
	}
	printf("\n");
		
	printf("\tParameters:  ");
	for(int i=0; i<params.size(); i++){
		if(i>0){
			printf(", ");
		}
		printf("%s",params[i].c_str());
	}
	printf("\n");
	
	//get rid of it (optional)
	if(demo){
		mModelFactory->deleteModel(m);
	}
}

//----------------------------------------------------------------------------------------------------------------------

#pragma mark Eigen matrix IO functions

void SaveMatrix(Eigen::MatrixXd matrix, std::string filename)
{
    std::ofstream f(filename.c_str());
	
	if(!f.is_open()){
		return;
	}
	
    int nrows = matrix.rows();
	int ncols = matrix.cols();
	f<<nrows<<"\n";
	f<<ncols<<"\n";
	for(int r=0; r<nrows; r++){
		for(int c=0; c<ncols; c++){
			f<<matrix(r,c)<<"\t";
		}
		f<<"\n";
	}
	f.close();
}

//---------------------------------------------------------------------------------------------------

Eigen::MatrixXd LoadMatrix(std::string filename)
{
	int nrows=0;
	int ncols=0;
	double value;
	Eigen::MatrixXd matrix;
	std::ifstream f(filename.c_str());
	if(!f.is_open()){
		return matrix;
	}
	f>>nrows;
	f>>ncols;
	matrix = Eigen::MatrixXd::Zero(nrows, ncols);;
	for(int r=0; r<nrows; r++){
		for(int c=0; c<ncols; c++){
			f>>value;
			matrix(r,c)=value;
		}
	}
	f.close();
	return matrix;
}	

