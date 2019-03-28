#pragma once
#include "basic.h"
#include <gurobi_c++.h>
#include "InputParam.h"

class SubSubLPProblem
{
private:
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar** f_ij_d;
	GRBVar* g_od;

	vector<GRBConstr> pi_ij;
	vector<GRBConstr> rau_ij_d;
	vector<GRBConstr> tau_od;

	int sind;
public:
	SubSubLPProblem(int scenarioIndex);

	SubSubDualSolution getOptimalSolution(SubMIPSolution &sms);
	SubSubDualSolution getOptimalSolution(SubMIPSolution_continuous &sms);

	~SubSubLPProblem()
	{
		InputParam *paramPointer = InputParam::getParam();
		InputParam param = *paramPointer;

		vector<GRBConstr>().swap(pi_ij);
		vector<GRBConstr>().swap(rau_ij_d);
		vector<GRBConstr>().swap(tau_od);

		for (int i = 0; i < param.legNum; i++)
			delete[] f_ij_d[i];
		delete[] g_od;
	}
};