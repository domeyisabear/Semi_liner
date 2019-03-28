#pragma once
#include "MasterProblem.h"
#include "SubLPProblem.h"
#include "SubMIPProblem.h"

class MP_benders :
	public MasterProblem
{
private:
	GRBEnv env =  GRBEnv();

	GRBModel model = GRBModel(env);

	GRBVar1d a_ij;
	GRBVar1d Z_w;

	vector<SubLPProblem*> slp_vector;
	vector<SubMIPProblem*> smp_vector;

public:
	MP_benders();
	MasterSolution getSolution() override;
	~MP_benders()
	{
		int len = slp_vector.size();
		for (int i = 0; i < len; i++)
			delete slp_vector[i];

		len = smp_vector.size();
		for (int i = 0; i < len; i++)
			delete smp_vector[i];
	}
};
