#pragma once
#include "basic.h"
#include <vector>

GRBVar1d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d);

GRBVar2d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1,int d2);


GRBVar3d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1, int d2, int d3);

vector<double> getOptimalDualValue(vector<GRBConstr>& constrs);
double cal_val(vector<double> &dual_val, vector<GRBConstr> &constr); // dual val * rhs

vector<int> getOptimalValue(GRBVar1d& var);
vector<vector<int>> getOptimalValue(GRBVar2d& var);