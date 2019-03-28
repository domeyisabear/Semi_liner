#pragma once
#include "basic.h"

class SubMIPProblem
{
public:
	virtual SubMIPSolution getOptimalSolution(MasterSolution &ms) = 0;
	virtual double getLowerBound() { return -MAX_NUM; };
	SubMIPProblem() {};
	virtual ~SubMIPProblem(){};
};

