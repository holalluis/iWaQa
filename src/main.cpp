/*
 *  main.cpp
 *  Main file for server (command interpreter)
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/MAIN
 *
 */ 
 
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "model.h"
#include "solver.h"
#include "evaluator.h"
#include "datatable.h"
#include "setup.h"
#include "server.h"
#include "sampleutils.h"

//############################################################################################################

iWQModelLayout * setup;		

//-----------------------------------------------------------------------------------------------------

void printHelp(std::string topics)
{
	std::transform(topics.begin(), topics.end(), topics.begin(), ::toupper);

	if(topics.size()==0){
		//prints general help and available commands
		printf("\nUsage:\n");
		printf("     - to run as a server, specify a valid layout filename.\n");
		printf("     - to execute a single command, supply a layout and the command with its parameters.\n");

		printf("\nAvailable commands:\n\n");
	}
	
	std::string act_cmd="";
	bool found=false;
	
	//RUN
	act_cmd="RUN";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("RUN - Run the model.\n");	
		printf("            Parameters:\n");
		printf("           (1) [parfile] initial parameter file\n");
		printf("               (optional, default=layout parameters)\n");
		printf("            2  [output_filename] output filename for results\n");
		printf("\n");
		found=true;
	}
	
	//RUN_UNCSIM
	act_cmd="RUN_UNCSIM";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("RUN_UNCSIM - Run the model with UNCSIM I/O format.\n");	
		printf("            Parameters:\n");
		printf("            1  [parfile] initial parameter file in UNCSIM format\n");
		printf("            2  [output_filename] output filename for results\n");
		printf("               (as UNCSIM model layout)\n");
		printf("\n");
		found=true;
	}
	
	//CAL
	act_cmd="CAL";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("CAL - Optimize parameters against the evaluation function.\n");	
		printf("            Parameters:\n");
		printf("           (1) [parfile] initial parameter file\n");
		printf("               (optional, default=layout parameters)\n");
		printf("            2  [output_filename] output filename for result parameters\n");
		printf("\n");
		found=true;
	}
	
	//EVAL
	act_cmd="EVAL";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("EVAL - Evaluate parameter set.\n");	
		printf("            Parameter:\n");
		printf("           (1) [parfile] parameter file (optional, default=layout parameters)\n");
		printf("\n");
		found=true;
	}
	
	//SENS_LOC
	act_cmd="SENS_LOC";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("SENS_LOC - Local sensitivity analysis.\n");	
		printf("            Parameters:\n");
		printf("            1  [target] name of the target variable\n");
		printf("           (2) [factor] perturbation factor (optional, default=0.1)\n");
		printf("            3  [output_filename] output filename for results\n");
		printf("\n");
		found=true;
	}
	
	//SENS_REG
	act_cmd="SENS_REG";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("SENS_REG - Variance-based regional sensitivity analysis.\n");	
		printf("            Parameters:\n");
		printf("            1  [target] name of the target variable\n");
		printf("            2  [factor] perturbation factor\n");
		printf("            3  [output_filename] output filename for results\n");
		printf("           (4) [numrounds] number of simulations (optional, default=500)\n");
		printf("\n");
		found=true;
	}
	
	//CONF_UNCSIM
	act_cmd="CONF_UNCSIM";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("CONF_UNCSIM - Translate layout to UNCSIM configuration.\n");	
		printf("            Parameter:\n");
		printf("            1  [output_dir] directory name for output files\n");
		printf("\n");
		found=true;
	}
	
	//INFO
	act_cmd="INFO";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("INFO - Get model information.\n");	
		printf("            Parameter:\n");
		printf("            1  [name] Identifier of a model instance or a type\n");
		printf("\n");
		found=true;
	}
	
	//MCMC
	act_cmd="MCMC";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("MCMC - Run Markov chain Monte Carlo sampling.\n");	
		printf("            Parameters:\n");
		printf("            1  [output_filename] output file for sample\n");
		printf("            2  [totallength] total number of rounds\n");
		printf("            3  [burninlength] number of rounds for burn-in\n");
		printf("           (4) [parameter_filename] starting parameter values (optional)\n");
		printf("           (5) [1/0] load proposal matrix? (optional)\n");
		printf("\n");
		found=true;
	}
	
	//MCMC_HAARIO
	act_cmd="MCMC_HAARIO";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("MCMC - Run Markov chain Monte Carlo sampling (Haario\'s method).\n");	
		printf("            Parameters:\n");
		printf("            1  [output_filename] output file for sample\n");
		printf("            2  [totallength] total number of rounds\n");
		printf("            3  [burninlength] number of rounds for burn-in\n");
		printf("           (4) [parameter_filename] starting parameter values (optional)\n");
		printf("\n");
		found=true;
	}
	
	//SAMPLE_HIST
	act_cmd="SAMPLE_HIST";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("SAMPLE_HIST - Create parameter distributions from MCMC samples.\n");	
		printf("            Parameters:\n");
		printf("            1  [sample_filename] MCMC sample file.\n");
		printf("            2  [output_filename] output file\n");
		printf("            3  [burninlength] number of rounds for burn-in\n");
		printf("\n");
		found=true;
	}
	
	//SAMPLE_CI
	act_cmd="RUN_SAMPLE";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("RUN_SAMPLE - Run on a parameter sample.\n");	
		printf("            Parameters:\n");
		printf("            1  [sample_filename] MCMC sample file.\n");
		printf("            2  [output_filename] output file\n");
		printf("\n");
		found=true;
	}
	
	//HELP
	act_cmd="HELP";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("HELP - Prints the above information.\n");	
		printf("            Parameter:\n");
		printf("           (1) [search_string] A string to find in commands\n");
		printf("               (optional, default prints everything).\n");
		printf("\n");
		found=true;
	}
	
	//NEW: SEQUENTIAL CALIBRATION & INPUT ADJUSTMENT COMMANDS
	act_cmd="SEQ_CAL";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("SEQ_CAL - Sequential parameter calibration procedure.\n");	
		printf("            Parameter:\n");
		printf("            1  [event_flag] Data field to separate events\n");
		printf("            2  [output_filename] output file\n");
		printf("\n");
		found=true;
	}
	
	act_cmd="SEQ_INP";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("SEQ_INP - Sequential input adjustment procedure.\n");	
		printf("            Parameter:\n");
		printf("            1  [event_flag] Data field to separate events\n");
		printf("            2  [input] Data field to calibrate\n");
		printf("            3  [output_filename] output file\n");
		printf("           (4) [priorname] prior for input (optional)\n");
		printf("\n");
		found=true;
	}
	//END NEW
	
	//DO_SERIES
	act_cmd="DO_SERIES";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("DO_SERIES - Generate model and error series from an MCMC sample.\n");	
		printf("            Parameter:\n");
		printf("           	1  [sample_filename] MCMC sample file.\n");
		printf("            2  Number of rounds to simulate from the sample.\n");
		printf("\n");
		found=true;
	}
	
	//DO_PRED_SERIES
	act_cmd="DO_PRED_SERIES";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("DO_PRED_SERIES - Generate model and error series from an MCMC sample (predictive mode only).\n");	
		printf("            Parameter:\n");
		printf("           	1  [sample_filename] MCMC sample file.\n");
		printf("            2  Number of rounds to simulate from the sample.\n");
		printf("\n");
		found=true;
	}
	
	//DO_BEST_SERIES
	act_cmd="DO_BEST_SERIES";
	if(topics.size()==0 || act_cmd.find(topics)!=std::string::npos){
		printf("DO_BEST_SERIES - Generate model and error series for a ML parameter file.\n");	
		printf("            Parameter:\n");
		printf("           	1  [parameter_filename] parameter file.\n");
		printf("\n");
		found=true;
	}
	
	if(!found && topics.size()){
		printf("No command found with \"%s\".\n",topics.c_str());
	}
	return;
}

//-----------------------------------------------------------------------------------------------------

std::string processcmd(std::string command)
{
	std::string answer="@I don't understand your command (\""+command+"\")\n";
	if(command.size()<2){
		//too short
		return answer;
	}
	else{
		if(command[0]!='@'){
			return answer;
		}
		if(command[command.size()-1]!='\n'){
			return answer;
		}
		command=command.substr(1,command.size()-2);
	}
	std::vector<std::string> tokens;
	Tokenize(command,tokens,"|");
	std::string pricommand="";
	
	if(tokens.size()>0){
		pricommand=tokens[0];
	}
	if(pricommand.compare("RUN")==0 && (tokens.size()==2 || tokens.size()==3)){
		if(setup->validity()<IWQ_VALID_FOR_RUN){
			answer="@Model layout is not valid for RUN.\n";
		}
		else{
			std::string parfilename;
			std::string outfilename;
			if(tokens.size()==2){
				outfilename=tokens[1];
			}
			if(tokens.size()==3){
				parfilename=tokens[1];
				outfilename=tokens[2];
				setup->loadParameters(parfilename);
			}
			setup->run();
			setup->saveResults(outfilename);
			return "@RUN completed.\n";
		}
	}
	if(pricommand.compare("RUN_UNCSIM")==0 && tokens.size()==3){
		if(setup->validity()<IWQ_VALID_FOR_RUN){
			answer="@Model layout is not valid for RUN.\n";
		}
		else{
			std::string parfilename;
			std::string outfilename;
			parfilename=tokens[1];
			outfilename=tokens[2];
			setup->loadParameters(parfilename, true);
			setup->run();
			setup->saveResultsUNCSIM(outfilename);
			return "@RUN_UNCSIM completed.\n";
		}
	}
	if(pricommand.compare("CAL")==0 && (tokens.size()==2 || tokens.size()==3)){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for CAL.\n";
		}
		else{
			//get filenames
			std::string paroutfilename;
			if(tokens.size()==3){
				//full config
				std::string parinfilename=tokens[1];
				paroutfilename=tokens[2];
				setup->loadParameters(parinfilename);
			}
			else{
				//just output filename given
				paroutfilename=tokens[1];
			}
			
			//do it
			setup->calibrate();
			setup->saveParameters(paroutfilename);
			return "@CAL completed.\n";
		}
	}
	if(pricommand.compare("EVAL")==0 && (tokens.size()==2 || tokens.size()==1)){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for EVAL.\n";
		}
		else{
			if(tokens.size()==2){
				std::string parfilename;
				parfilename=tokens[1];
				setup->loadParameters(parfilename);
			}
			double evalresult=setup->evaluate();
			std::stringstream s;
			s<<"@EVAL returned "<<evalresult<<"\n";
			return s.str();
		}
	}
	if(pricommand.compare("SENS_LOC")==0 && (tokens.size()==3 || tokens.size()==4)){
		if(setup->validity()<IWQ_VALID_FOR_RUN){
			answer="@Model layout is not valid for SENS_LOC.\n";
		}
		else{
			std::string target=tokens[1];
			double factor=0.1;
			std::string filename;
			if(tokens.size()==4){
				char * endptr;
				factor=strtod(tokens[2].c_str(), &endptr);
				if(*endptr || factor<=0.0){
					return "@Perturbation factor is not a valid number.\n";
				}
				filename=tokens[3];
			}
			else{
				filename=tokens[2];
			}
			setup->localSensitivityAnalysis(factor,target,filename);
			
			return "@SENS_LOC completed.\n";
		}
	}
	if(pricommand.compare("SENS_REG")==0 && (tokens.size()==4 || tokens.size()==5)){
		if(setup->validity()<IWQ_VALID_FOR_RUN){
			answer="@Model layout is not valid for SENS_REG.\n";
		}
		else{
			std::string target=tokens[1];
			double factor=0.1;
			std::string filename;
			char * endptr;
			factor=strtod(tokens[2].c_str(), &endptr);
			if(*endptr || factor<=0.0){
				return "@Perturbation factor is not a valid number.\n";
			}
			filename=tokens[3];
			int numtrials=500;
			if(tokens.size()==5){
				numtrials=strtol(tokens[4].c_str(), &endptr, 10);
				if(*endptr || numtrials<3){
					return "@Sample size is not a valid number or less than 3.\n";
				}
			}
			setup->regionalSensitivityAnalysis(factor,target,filename,numtrials);
			
			return "@SENS_REG completed.\n";
		}
	}
	if(pricommand.compare("CONF_UNCSIM")==0 && tokens.size()==2){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for CONF_UNCSIM.\n";
		}
		else{
			std::string target=tokens[1];
			setup->furnishUNCSIM(target);
			return "@CONF_UNCSIM completed.\n";
		}
	}
	if(pricommand.compare("INFO")==0 && tokens.size()==2){
		std::string name=tokens[1];
		setup->printModelInfo(name);
		return "@INFO completed.\n";
	}
	if(pricommand.compare("MCMC")==0 && (tokens.size()==4 || tokens.size()==5 || tokens.size()==6)){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for MCMC.\n";
		}
		else{
			std::string filename=tokens[1];
			int numrounds=5500;
			int burnin=500;
			int loadprop = 0;
					
			if(tokens.size()==5 || tokens.size()==6){
				std::string parfilename = tokens[4];
				//set these values
				setup->loadParameters(parfilename);
			}
			
			char * endptr;
			if(tokens.size()==6){
				loadprop = strtol(tokens[5].c_str(), &endptr, 0);
			}
			
			numrounds=strtol(tokens[2].c_str(), &endptr, 10);
			if(*endptr || numrounds<=0){
				return "@Iteration count is not a valid number.\n";
			}
			
			burnin=strtol(tokens[3].c_str(), &endptr, 10);
			if(*endptr || burnin<0){
				return "@Burn-in length is not a valid number.\n";
			}
			
			setup->MCMC(numrounds, burnin, filename, (bool)loadprop);
			
			return "@MCMC completed.\n";
		}
	}
	if(pricommand.compare("MCMC_HAARIO")==0 && (tokens.size()==4 || tokens.size()==5)){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for MCMC_HAARIO.\n";
		}
		else{
			std::string filename=tokens[1];
			int numrounds=5500;
			int burnin=500;
			
			if(tokens.size()==5){
				std::string parfilename = tokens[4];
				//set these values
				setup->loadParameters(parfilename);
			}
			
			char * endptr;
			numrounds=strtol(tokens[2].c_str(), &endptr, 10);
			if(*endptr || numrounds<=0){
				return "@Iteration count is not a valid number.\n";
			}
			
			burnin=strtol(tokens[3].c_str(), &endptr, 10);
			if(*endptr || burnin<0){
				return "@Burn-in length is not a valid number.\n";
			}
			
			setup->MCMC_Haario(numrounds, burnin, filename);
			
			return "@MCMC completed.\n";
		}
	}
	if(pricommand.compare("SAMPLE_HIST")==0 && tokens.size()==4){
		std::string filename=tokens[1];
		std::string ofilename=tokens[2];
		int burnin=0;
		
		char * endptr;
		burnin=strtol(tokens[3].c_str(), &endptr, 10);
		if(*endptr || burnin<0){
			return "@Burn-in length is not a valid number.\n";
		}
			
		iWQEstimateDistributions(filename, setup->parameters()->namesForPlainValues(), ofilename, burnin);
			
		return "@SAMPLE_HIST completed.\n";
	}
	if(pricommand.compare("RUN_SAMPLE")==0 && tokens.size()==3){
		std::string filename=tokens[1];
		std::string ofilename=tokens[2];
				
		setup->runOnSample(filename, ofilename);	
			
		return "@RUN_SAMPLE completed.\n";
	}
	if(pricommand.compare("DO_SERIES")==0 && tokens.size()==3){
		std::string filename=tokens[1];
		char * endptr;
		int rounds=strtol(tokens[2].c_str(), &endptr, 10);
		if(*endptr || rounds<=0){
			return "@Desired rounds is not a valid number.\n";
		}
				
		setup->runStandardSeriesOnSample(filename, rounds, false, false);
			
		return "@DO_SERIES completed.\n";
	}	
	if(pricommand.compare("DO_PRED_SERIES")==0 && tokens.size()==3){
		std::string filename=tokens[1];
		char * endptr;
		int rounds=strtol(tokens[2].c_str(), &endptr, 10);
		if(*endptr || rounds<=0){
			return "@Desired rounds is not a valid number.\n";
		}
				
		setup->runStandardSeriesOnSample(filename, rounds, true, false);
			
		return "@DO_PRED_SERIES completed.\n";
	}
	if(pricommand.compare("DO_BEST_SERIES")==0 && tokens.size()==2){
		std::string filename=tokens[1];
		setup->createBestSeries(filename);
			
		return "@DO_BEST_SERIES completed.\n";
	}
	if(pricommand.compare("HELP")==0){
		std::string searchstr="";
		if(tokens.size()>1){
			searchstr=tokens[1];
		}
		printHelp(searchstr);		
		return "@HELP completed.\n";
	}
	//NEW: SPECIFIC EVENT-DEPENDENT OPTIMISATION PROCEDURES
	if(pricommand.compare("SEQ_CAL")==0 && tokens.size()==3){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for SEQ_CAL.\n";
		}
		else{
			std::string eventflag=tokens[1];
			std::string target=tokens[2];
			setup->evaluator()->sequentialCalibrateParameters(eventflag);
			setup->dataTable()->writeToFile(target);
			return "@SEQ_CAL completed.\n";
		}
	}
	
	if(pricommand.compare("SEQ_INP")==0 && (tokens.size()==4 || tokens.size()==5)){
		if(setup->validity()<IWQ_VALID_FOR_CALIBRATE){
			answer="@Model layout is not valid for SEQ_INP.\n";
		}
		else{
			std::string eventflag=tokens[1];
			std::string inpname=tokens[2];
			std::string target=tokens[3];
			iWQRandomGenerator * prior=NULL;
			if(tokens.size()==5){
				prior=setup->distributionForName(tokens[4]);
			}
			setup->evaluator()->sequentialCalibrateInputs(eventflag, inpname, prior);
			setup->dataTable()->writeToFile(target);
			return "@SEQ_INP completed.\n";
		}
	}
	//END NEW
	
	return answer;
}

//############################################################################################################

int main (int argc, char * const argv[])
{
	if(argc<2){
		printf("[Error]: Specify the layout filename or type \"HELP\".\n");
		
		//TODO: make exceptions with those commands which do not require a layout (HELP, INFO*, SAMPLE_HIST)
		
		return 1;
	}
	
	std::string argv1=argv[1];
	if(argv1.compare("HELP")==0){
		std::string searchstr="";
		if(argc>2){
			searchstr=argv[2];
		}
		printHelp(searchstr);
		return 0;
	}
	
	setup=new iWQModelLayout(argv1);
	
	if(setup->validity()<IWQ_VALID_FOR_RUN){
		printf("Model setup is invalid to run.\n");
		return 1;
	}
	
	//export layout
	setup->saveLayoutGraph("model_layout.dot");
		
	// Run-on-demand
	if(argc==2){
		iWQServer::run(5555,processcmd);
	}
	else{
		printf("*** Offline mode ***\n");
		std::string command="@";
		for(int i=2; i<argc; i++){
			if(i>2){
				command+="|";
			}
			command+=argv[i];
		}
		command+="\n";
		std::string result=processcmd(command);
		result=result.substr(1,result.size()-1);
		printf("Result: %s\n",result.c_str());
	}
	
	return 0;
}

//############################################################################################################

