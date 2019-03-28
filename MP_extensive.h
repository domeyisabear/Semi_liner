#pragma once
#include "basic.h"
#include "MasterProblem.h"
#include "InputParam.h"

class MP_extensive :
	public MasterProblem
{
private:
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar1d a_ij;
	GRBVar2d b_w_ij;
	GRBVar3d x_w_ij_s;
	GRBVar2d y_w_i;
	GRBVar2d f_w_ij;
	GRBVar2d fv_w_ij;
	GRBVar2d g_w_h;
	GRBVar2d gama_w_h;
	GRBVar2d t_w_i;

	vector<GRBLinExpr> expr_v;

	bool isSolved;

public:
	MasterSolution getSolution() override;
	MP_extensive();
	~MP_extensive()
	{
	}
};
