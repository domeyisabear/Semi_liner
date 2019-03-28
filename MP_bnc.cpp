#include "MP_bnc.h"
#include "SMP_extensive.h"
#include "tools.h"
#include "MPCallback.h"

MP_bnc::MP_bnc()
{
	InputParam *paramPointer = InputParam::getParam();
	InputParam param = *paramPointer;

	a_ij = addVars(model, 0, 1, GRB_BINARY, "a", param.fixPrioList.size());
	Z_w = addVars(model, -MAX_NUM, GRB_INFINITY, GRB_CONTINUOUS, "Z", param.scenarioNum);

	// constraints
	// constraint (1)
	{
		set<set<int>> S;
		for (int i = 0; i < param.fixPrioList.size(); i++)
		{
			int id_i = param.fixPrioList[i].first;
			int id_j = param.fixPrioList[i].second;
			int prio2 = param.getFixPrioIndex(id_j, id_i);

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
		for (int i = 0; i < param.fixPortIndices.size(); i++)
		{
			for (int j = 0; j < param.fixPortIndices.size(); j++)
			{
				if (i == j)
					continue;
				for (int k = 0; k < param.fixPortIndices.size(); k++)
				{
					if ((i == k) || (j == k))
						continue;
					int ind_i = param.fixPortIndices[i];
					int ind_j = param.fixPortIndices[j];
					int ind_k = param.fixPortIndices[k];

					int port1 = param.portList[ind_i].portID;
					int port2 = param.portList[ind_j].portID;
					int port3 = param.portList[ind_k].portID;

					int index1 = param.getFixPrioIndex(port1, port2);
					int index2 = param.getFixPrioIndex(port2, port3);
					int index3 = param.getFixPrioIndex(port3, port1);

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
	for (int w = 0; w < param.scenarioNum; w++)
	{
		obj += param.scenarioList[w].probability*Z_w[w];
	}

	model.setObjective(obj, GRB_MINIMIZE);

	//model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	slp_vector = vector<SubLPProblem*>(param.scenarioNum);
	smp_vector = vector<SubMIPProblem*>(param.scenarioNum);

	for (int w = 0; w < param.scenarioNum; w++)
	{
		slp_vector[w] = new SubLPProblem(w);
		smp_vector[w] = new SMP_extensive(w);
		Z_w[w].set(GRB_DoubleAttr_LB,smp_vector[w]->getLowerBound());
	}

	mpcb = new MPCallback(a_ij, Z_w, slp_vector, smp_vector);

	model.setCallback(mpcb);

	model.set(GRB_IntParam_LazyConstraints, 1);
	//model.getEnv().set(GRB_IntParam_DualReductions, 0);
}

MasterSolution MP_bnc::getSolution()
{
	model.optimize();

	MasterSolution ms;

	InputParam *paramPointer = InputParam::getParam();
	InputParam param = *paramPointer;

	double *val_a_ij = model.get(GRB_DoubleAttr_X, a_ij, param.fixPrioList.size());
	ms.a_ij = vector<int>(param.fixPrioList.size());
	for (int i = 0; i < param.fixPrioList.size(); i++)
		ms.a_ij[i] = round(val_a_ij[i]);
	free(val_a_ij);

	double *val_Z_w = model.get(GRB_DoubleAttr_X, Z_w, param.scenarioNum);
	ms.Z_w = vector<double>(param.scenarioNum);
	for (int i = 0; i < param.scenarioNum; i++)
		ms.Z_w[i] = val_Z_w[i];
	free(val_Z_w);
	ms.objVal = model.get(GRB_DoubleAttr_ObjVal);

	return ms;
}
