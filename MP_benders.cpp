#include "MP_benders.h"
#include "tools.h"
#include "SMP_extensive.h"
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

MP_benders::MP_benders()
{
	a_ij = addVars(model, 0, 1, GRB_BINARY, "a", param->fixPrioList.size());
	Z_w = addVars(model, -MAX_NUM, GRB_INFINITY, GRB_CONTINUOUS, "Z", param->scenarioNum);

	// constraints
	// constraint (1)
	{
		set<set<int>> S;
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			int id_i = param->fixPrioList[i].first;
			int id_j = param->fixPrioList[i].second;
			int prio2 = param->getFixPrioIndex(id_j, id_i);

			set<int> s;
			s.insert(i);
			s.insert(prio2);

			int len1 = S.size();
			S.insert(s);
			int len2 = S.size();

			if (len2>len1)
				model.addConstr(a_ij[i] + a_ij[prio2] == 1);
		}
	}

	// constraint (2)
	{
		set<set<int>> S;
		for (int i = 0; i < param->fixPortIndices.size(); i++)
		{
			for (int j = 0; j < param->fixPortIndices.size(); j++)
			{
				if (i == j)
					continue;
				for (int k = 0; k < param->fixPortIndices.size(); k++)
				{
					if ((i == k) || (j == k))
						continue;
					int ind_i = param->fixPortIndices[i];
					int ind_j = param->fixPortIndices[j];
					int ind_k = param->fixPortIndices[k];

					int port1 = param->portList[ind_i].portID;
					int port2 = param->portList[ind_j].portID;
					int port3 = param->portList[ind_k].portID;

					int index1 = param->getFixPrioIndex(port1, port2);
					int index2 = param->getFixPrioIndex(port2, port3);
					int index3 = param->getFixPrioIndex(port3, port1);

					set<int> s;
					s.insert(index1);
					s.insert(index2);
					s.insert(index3);

					int len1 = S.size();
					S.insert(s);
					int len2 = S.size();

					if (len2>len1)
						model.addConstr(a_ij[index1] + a_ij[index2] + a_ij[index3] <= 2);
				}
			}
		}
	}

	// objective
	GRBLinExpr obj = 0;
	for (int w = 0; w < param->scenarioNum; w++)
	{
		obj += param->scenarioList[w].probability*Z_w[w];
	}

	model.setObjective(obj, GRB_MINIMIZE);

	model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	int thread_num = std::thread::hardware_concurrency();
	if (thread_num)
		model.set(GRB_IntParam_Threads, thread_num);

	slp_vector = vector<SubLPProblem*>(param->scenarioNum);
	smp_vector = vector<SubMIPProblem*>(param->scenarioNum);

	std::vector< std::future<void> > results;
	
	for (int w = 0; w < param->scenarioNum; w++)
	{
		results.push_back(
			param->executor->commit([](vector<SubLPProblem*> *slbp_v, vector<SubMIPProblem*> *smp_v, GRBVar1d Z, int i)
		{
			(*slbp_v)[i] = new SubLPProblem(i); (*smp_v)[i] = new SMP_extensive(i); Z[i].set(GRB_DoubleAttr_LB, (*smp_v)[i]->getLowerBound());
		}, &slp_vector, &smp_vector, Z_w, w)
			);
	}
	for (int w = 0; w < param->scenarioNum; w++)
	{
		results[w].get();
	}
}

MasterSolution MP_benders::getSolution()
{
	clock_t start = clock();

	double UB = MAX_NUM;
	double LB = -MAX_NUM;

	ofstream outputfile(param->getWorkDir() + "UB_LB_benders.csv");

	model.optimize();

	LB = model.get(GRB_DoubleAttr_ObjVal);
	std::cout << endl;
	cout << "----------------------------------------" << endl;
	cout << setprecision(12) << "UB : " << UB << " LB : " << LB << ", 0" << endl;
	cout << "diff : " << (UB - LB) / fabs(LB + 1e-9) * 100 << "%" << endl;
	outputfile << setprecision(16) << UB << "," << LB << "," << 0 << "\n";

	int iter_num = 0;

	while (UB - LB>1e-3)
	{
		MasterSolution ms;

		ms.a_ij = vector<int>(param->fixPrioList.size());
		for (int i = 0; i < param->fixPrioList.size(); i++)
		{
			double val = a_ij[i].get(GRB_DoubleAttr_X);
			ms.a_ij[i] = Round(val);
		}

		ms.Z_w = vector<double>(param->scenarioNum);
		for (int i = 0; i < param->scenarioNum; i++)
			ms.Z_w[i] = Z_w[i].get(GRB_DoubleAttr_X);
		
		ms.objVal = model.get(GRB_DoubleAttr_ObjVal);

		bool violation = false;

		std::vector< std::future<SubDualSolution> > results_slp;
		for (int w = 0; w < param->scenarioNum; w++)
		{
			results_slp.push_back(param->executor->commit([](SubLPProblem* slp, MasterSolution &ms_) {return slp->getOptimalSolution(ms_); }, slp_vector[w], ms));
		}

		for (int w = 0; w < param->scenarioNum; w++)
		{
			SubDualSolution sds = results_slp[w].get();
			if (sds.objVal - ms.Z_w[w] >= 1) // add lp cut
			{
				violation = true;
				int iter = 0;
				GRBLinExpr expr = 0;
				for (int i = 0; i < param->fixPrioList.size(); i++)
				{
					expr += a_ij[i] * sds.xi_ij[i];
					iter++;
				}
				expr += sds.other_const;
				GRBConstr constr = model.addConstr(Z_w[w] - expr >= 0);
				constr.set(GRB_IntAttr_Lazy, 1);
			}
			sds.end();
		}

		if (violation)
		{
			model.optimize();
			LB = model.get(GRB_DoubleAttr_ObjVal);
			ms.end();

			// check time limit
			double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
			if (time_elps>timeLimit)
			{
				std::cout << "UB: " << UB << ", LB:" << LB << ", 0" << ", time elapsed: " << time_elps << endl;
				outputfile << setprecision(16) << UB << "," << LB << "," << 0 << "\n";
				break;
			}

			cout << "-------------------------------------" << endl;
			std::cout << "UB: " << UB << ", LB:" << LB << ", 0" << ", time elapsed: " << time_elps << endl;
			cout << "diff : " << (UB - LB) / fabs(LB + 1e-9) * 100 << "%" << endl;
			outputfile << setprecision(16) << UB << "," << LB << "," << 0 << "\n";
			continue;
		}

		GRBLinExpr expr_1 = 0;
		GRBLinExpr expr_0 = 0;
		int num = 0;

		for (int i = 0; i < ms.a_ij.size(); i++)
		{
			if (ms.a_ij[i] == 1)
			{
				expr_1 += a_ij[i];
				num++;
			}
			else
				expr_0 += a_ij[i];
		}

		double t_UB = 0;

		std::vector< std::future<SubMIPSolution> > results_smp;
		for (int w = 0; w < param->scenarioNum; w++)
		{
			results_smp.push_back(param->executor->commit([](SubMIPProblem* smp, MasterSolution &ms_) {return smp->getOptimalSolution(ms_); }, smp_vector[w], ms));
		}

		for (int w = 0; w < param->scenarioNum; w++)
		{
			//cout << w << endl;
			SubMIPSolution sms = results_smp[w].get();

			if(sms.objVal==MAX_NUM)
			{
				t_UB = MAX_NUM;
				GRBConstr constr = model.addConstr(num - expr_1 + expr_0 >= 2);
				constr.set(GRB_IntAttr_Lazy, 1);

				for (w++; w < param->scenarioNum; w++)
					results_smp[w].wait();

				break;
			}

			t_UB += param->scenarioList[w].probability*sms.objVal;
			//cout << sms.objVal << endl;
			if (sms.objVal - ms.Z_w[w] >= 1e-6) // add mip cut
			{
				double L = smp_vector[w]->getLowerBound();
				GRBConstr constr = model.addConstr(Z_w[w] - ((sms.objVal - L)*(expr_1 - expr_0 - num + 1) + L) >= 0);
				constr.set(GRB_IntAttr_Lazy, 1);
			}
			sms.end();
		}

		if (UB > t_UB)
			UB = t_UB;

		model.optimize();
		LB = model.get(GRB_DoubleAttr_ObjVal);
		ms.end();

		// check time limit
		double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
		if (time_elps>timeLimit)
		{
			std::cout << "UB: " << UB << ", LB:" << LB << ", 1" << ", time elapsed: " << time_elps << endl;
			outputfile << setprecision(16) << UB << "," << LB << "," << 1 << "\n";
			break;
		}

		// check iteration limit
		iter_num++;
		if (iterationLim > 0 && iter_num >= iterationLim)
		{
			break;
		}

		cout << endl;
		cout << "-------------------------------------" << endl;
		std::cout << "UB: " << UB << ", LB:" << LB << ", 1" << ", time elapsed: " << time_elps << endl;
		cout << "diff : " << (UB - LB) / fabs(LB + 1e-9) * 100 << "%" << endl;
		outputfile << setprecision(16) << UB << "," << LB << "," << 1 << "\n";
	}

	outputfile.close();

	MasterSolution ms;

	ms.a_ij = vector<int>(param->fixPrioList.size());
	for (int i = 0; i < param->fixPrioList.size(); i++)
	{
		double val = a_ij[i].get(GRB_DoubleAttr_X);
		ms.a_ij[i] = Round(val);
	}

	ms.Z_w = vector<double>(param->scenarioNum);
	for (int i = 0; i < param->scenarioNum; i++)
		ms.Z_w[i] = Z_w[i].get(GRB_DoubleAttr_X);

	ms.objVal = model.get(GRB_DoubleAttr_ObjVal);

	return ms;
}
