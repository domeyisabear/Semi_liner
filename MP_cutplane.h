#pragma once
#include "MasterProblem.h"
#include "InputParam.h"
#include "SubMIPProblem.h"
#include "SubLBProblem.h"
#include "SubLBProblem_benders.h"

class MP_cutplane :
	public MasterProblem
{
private:
	
	vector<SubLBProblem*> slbp_vector;
	vector<SubMIPProblem*> smp_vector;

public:
	MP_cutplane();
	MasterSolution getSolution() override;
	~MP_cutplane()
	{
		for (int w = 0; w < param->scenarioNum;w++)
		{
			delete slbp_vector[w];
			delete smp_vector[w];
		}
	}
};

