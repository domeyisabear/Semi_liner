#pragma once
#include "MasterProblem.h"
#include "InputParam.h"
#include "SubMIPProblem.h"
#include "SubLBProblem.h"

class MP_cutplane_multi :
	public MasterProblem
{
private:

	vector<SubLBProblem*> slbp_vector;
	vector<SubMIPProblem*> smp_vector;

public:
	MP_cutplane_multi();
	MasterSolution getSolution() override;
	~MP_cutplane_multi()
	{
		for (int w = 0; w < param->scenarioNum; w++)
		{
			delete slbp_vector[w];
			delete smp_vector[w];
		}
	}
};

