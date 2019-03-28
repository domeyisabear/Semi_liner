#pragma once
#include "basic.h"
#include "gurobi_c++.h"
#include <ctime>

class SLBPCallback :
	public GRBCallback
{
private:
	double time_lim;
	clock_t start;
	double bestObjVal;

protected:
	void callback() override;
public:
	SLBPCallback(double timeLimit);
	void reset();
	~SLBPCallback(){}
};