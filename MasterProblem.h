#pragma once
#include "basic.h"

class MasterProblem
{
public:
	virtual MasterSolution getSolution()=0;
	virtual ~MasterProblem(){};
};