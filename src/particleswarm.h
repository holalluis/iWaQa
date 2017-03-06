/*
 *  particleswarm.h
 *  Particle-swarm optimiser
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/OPTIMISE
 *
 */ 

#include <vector>
#include <string>

#ifndef particleswarm_h
#define particleswarm_h

class iWQEvaluator;

typedef double (*iWQEvalFunction)(double *, int);

class iWQBounds
{
public:
	double min;
	double max;
};

iWQBounds iWQMakeBounds(double _min, double _max);

class iWQBoundsList : public std::vector<iWQBounds>
{
public:
	void add(iWQBounds b){ push_back(b); }
	void add(double min, double max);
};

//----------------------------------------------------------------------------

std::vector<double> iWQParticleSwarmOptimize(iWQEvaluator * evaluator, iWQBoundsList bounds, int populationsize, int maxiterations=2000, int idlerunlength=50);

#endif
