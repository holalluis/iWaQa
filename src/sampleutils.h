/*
 *  sampleutils.h
 *  Histogram creator from MCMC samples
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/UTILS
 *
 */ 

#include <string>
#include <vector>

#ifndef sampleutils_h
#define sampleutils_h
  
void iWQEstimateDistributions(std::string samplefilename, iWQStrings paramnames, std::string outputfilename, int burn_in_length);

#endif

