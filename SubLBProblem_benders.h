#pragma once
#include <gurobi_c++.h>
#include "MasterProblem.h"
#include "InputParam.h"
class SubLBProblem_benders
{
private:
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar1d b_ij;
	GRBVar2d x_ij_s;
	GRBVar1d y_i;
	GRBVar1d gama_h;
	GRBVar1d t_i;

	GRBVar B;

	int sind;

	class check_feasible
	{
	private:
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBVar1d a_ij;
	public:
		check_feasible();
		void removeSolution(MasterSolution&);
		void removeSolutions(vector<MasterSolution>&);
		bool check();
		~check_feasible() {};
	};

	check_feasible* feasible = new check_feasible();

public:
	SubLBProblem_benders(int scenarioIndex);
	MasterSolution getSolution();
	//vector<MasterSolution> getSolutions();
	void setDualPenalty(vector<double>& omega);
	void removeSolution(MasterSolution&);
	void removeSolutions(vector<MasterSolution>&);
	~SubLBProblem_benders()
	{
		delete feasible;
	}
};

