#pragma once
#include <gurobi_c++.h>
#include "basic.h"
#include "InputParam.h"

bool subTourSep(GRBVar* x_ij, double* val_x_ij, GRBVar* y_i, double* val_y_i, vector<GRBConstr>& cut);

class SubLPProblem
{
private:
	GRBEnv env =  GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar1d b_ij;
	GRBVar2d x_ij_s;
	GRBVar1d y_i;
	GRBVar1d f_ij;
	GRBVar1d fv_ij;
	GRBVar1d g_h;
	GRBVar1d gama_h;
	GRBVar1d t_i;

	vector<GRBConstr> a_constr;
	vector<GRBConstr> other_constr;

	int sind;

public:
	SubLPProblem(int scenarioIndex);
	SubDualSolution getOptimalSolution(MasterSolution &ms);
	SubDualSolution getOptimalSolution(MasterSolution_continuous &ms);
	~SubLPProblem()
	{};
};

