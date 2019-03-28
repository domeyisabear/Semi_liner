#pragma once
#include "basic.h"
#include "InputParam.h"
#include "tools.h"

class LinerProblem
{
private:
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar* a_ij;
	GRBVar* x_ij;
	GRBVar* y_i;
	GRBVar** f_w_ij;
	GRBVar** g_w_od;

	GRBLinExpr obj = 0;
	GRBLinExpr legCost = 0;
	GRBLinExpr portCost = 0;

	vector<GRBLinExpr> frtRev;
	vector<GRBLinExpr> frtAmt;

public:
	LinerProblem();
	void getSolution();
	~LinerProblem()
	{
		delete[] a_ij;
		delete[] x_ij;
		delete[] y_i;
		for (int i = 0; i < param->scenarioNum; i++)
		{
			delete[] g_w_od[i];
			delete[] f_w_ij[i];
		}
		delete[] f_w_ij;
		delete[] g_w_od;
	}
};

