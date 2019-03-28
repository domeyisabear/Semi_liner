#pragma once
#include "basic.h"
#include "tools.h"
#include <vector>

GRBVar1d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d)
{
	vector<GRBVar> vars(d);
	for (int i = 0; i < d;i++)
	{
		string str=v_name+"_"+to_string(i);
		vars[i] = model.addVar(lb, ub, 0, vtype, str);
	}
	return vars;
}

GRBVar2d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1, int d2)
{
	GRBVar2d vars(d1);
	for (int i = 0; i < d1;i++)
	{
		vars[i] = GRBVar1d(d2);
		for (int j = 0; j < d2;j++)
		{
			string str=v_name+"_"+to_string(i)+"_"+to_string(j);
			vars[i][j] = model.addVar(lb, ub, 0, vtype, str);
		}
	}
	return vars;
}

GRBVar3d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1, int d2, int d3)
{
	GRBVar3d vars(d1);
	for (int i = 0; i < d1;i++)
	{
		vars[i] = GRBVar2d(d2);
		for (int j = 0; j < d2;j++)
		{
			vars[i][j] = GRBVar1d(d3);
			for (int k = 0; k < d3;k++)
			{
				string str=v_name+"_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k);
				vars[i][j][k] = model.addVar(lb, ub, 0, vtype, str);
			}
		}
	}
	return vars;
}

vector<double> getOptimalDualValue(vector<GRBConstr>& constrs)
{
	int len = constrs.size();
	vector<double> re = vector<double>(len);
	for (int i = 0; i < len; i++)
		re[i] = constrs[i].get(GRB_DoubleAttr_Pi);
	return re;
}

double cal_val(vector<double> &dual_val, vector<GRBConstr> &constr)
{
	double tmp = 0;
	for (int i = 0; i < constr.size(); i++)
	{
		tmp += dual_val[i] * constr[i].get(GRB_DoubleAttr_RHS);
	}
	return tmp;
}

vector<int> getOptimalValue(GRBVar1d& var)
{
	int len = var.size();
	vector<int> val = vector<int>(len);
	for (int i = 0; i < len;i++)
	{
		double a= var[i].get(GRB_DoubleAttr_X);
		val[i] = Round(a);
	}

	return val;
}

vector<vector<int>> getOptimalValue(GRBVar2d& var)
{
	int len = var.size();
	vector<vector<int>> val = vector<vector<int>>(len);
	for (int i = 0; i < len; i++)
	{
		val[i] = vector<int>(var[i].size());
		for (int j = 0; j < var[i].size(); j++)
		{
			double a = var[i][j].get(GRB_DoubleAttr_X);
			val[i][j] = Round(a);
		}
	}

	return val;
}