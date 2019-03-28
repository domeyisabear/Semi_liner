#include "SubLBProblem_benders.h"
#include "tools.h"

SubLBProblem_benders::check_feasible::check_feasible()
{
	a_ij = addVars(model, 0, 1, GRB_BINARY, "a", param->fixPrioList.size());

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

			if (len2 > len1)
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

					if (len2 > len1)
						model.addConstr(a_ij[index1] + a_ij[index2] + a_ij[index3] <= 2);
				}
			}
		}
	}

	model.setObjective(GRBLinExpr(0), GRB_MINIMIZE);
	model.getEnv().set(GRB_IntParam_OutputFlag, 0);
}

void SubLBProblem_benders::check_feasible::removeSolution(MasterSolution& ms)
{
	GRBLinExpr expr = 0;
	for (int i = 0; i < param->fixPrioList.size(); i++)
	{
		if (ms.a_ij[i] == 1)
			expr += 1 - a_ij[i];
		else
			expr += a_ij[i];
	}
	model.addConstr(expr >= 2);////////////////////////////////////
}

void SubLBProblem_benders::check_feasible::removeSolutions(vector<MasterSolution>& vector_ms)
{
	for (int j = 0; j < vector_ms.size(); j++)
	{
		removeSolution(vector_ms[j]);
	}
}

bool SubLBProblem_benders::check_feasible::check()
{
	model.optimize();
	int result = model.get(GRB_IntAttr_Status);
	if (result == GRB_INFEASIBLE || result == GRB_INF_OR_UNBD)
		return false;
	else
		return true;
}

SubLBProblem_benders::SubLBProblem_benders(int scenarioIndex) :sind(scenarioIndex)
{
	vector<Demand> demandList = param->scenarioList[scenarioIndex].demandList;

	int demandNum = demandList.size();

	b_ij = addVars(model, 0, 1, GRB_BINARY, "b", param->prioList.size());
	x_ij_s = GRBVar2d(param->legNum);
	for (int a = 0; a < param->legNum; a++)
	{
		string str = "x_" + to_string(a);
		x_ij_s[a] = addVars(model, 0, 1, GRB_BINARY, str, param->legList[a].confList.size());
	}
	y_i = addVars(model, 0, 1, GRB_BINARY, "y", param->addPortIndices.size());
	gama_h = addVars(model, 0, 1, GRB_BINARY, "gama", demandNum);
	t_i = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "t", param->portNum);

	B = model.addVar(-MAX_NUM, GRB_INFINITY, 0, GRB_CONTINUOUS, "B");

	// constraint (3)
	{
		set<set<int>> S;
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;

			int prio2 = param->getPrioIndex(id_j, id_i);

			set<int> s;
			s.insert(i);
			s.insert(prio2);

			int len1 = S.size();
			S.insert(s);
			int len2 = S.size();

			if (len2 > len1)
				model.addConstr(b_ij[i] + b_ij[prio2] == 1);
		}
	}

	// constraint (4)
	{
		set<set<int>> S;
		for (int i = 0; i < param->portNum; i++)
		{
			if (i == param->baseStartIndex || i == param->baseEndIndex)
				continue;
			for (int j = 0; j < param->portNum; j++)
			{
				if (j == param->baseStartIndex || j == param->baseEndIndex)
					continue;
				if (i == j)
					continue;
				for (int k = 0; k < param->portNum; k++)
				{
					if (k == param->baseStartIndex || k == param->baseEndIndex)
						continue;
					if (j == k || k == i)
						continue;
					int port_i = param->portList[i].portID;
					int port_j = param->portList[j].portID;
					int port_k = param->portList[k].portID;
					int prio_ij = param->getPrioIndex(port_i, port_j);
					int prio_jk = param->getPrioIndex(port_j, port_k);
					int prio_ki = param->getPrioIndex(port_k, port_i);
					if (prio_ij == -1 || prio_jk == -1 || prio_ki == -1)
						continue;

					set<int> s;
					s.insert(prio_ij);
					s.insert(prio_jk);
					s.insert(prio_ki);

					int len1 = S.size();
					S.insert(s);
					int len2 = S.size();

					if (len2 > len1)
						model.addConstr(b_ij[prio_ij] + b_ij[prio_jk] + b_ij[prio_ki] <= 2);
				}
			}
		}
	}

	// constraint (8)
	{
		for (int a = 0; a < param->legNum; a++)
		{
			GRBLinExpr expr = 0;
			for (GRBVar& var : x_ij_s[a])
			{
				expr += var;
			}
			model.addConstr(expr <= 1);
		}
	}

	// constraint (9)
	{
		// compulsory port
		for (int i = 0; i < param->fixPortIndices.size(); i++)
		{
			GRBLinExpr in_expr = 0;
			GRBLinExpr out_expr = 0;
			int ind_i = param->fixPortIndices[i];
			int port_i = param->portList[ind_i].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				int out_index = param->getLegIndex(port_i, port_j);
				if (in_index != -1)
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr == 1);
			model.addConstr(in_expr == 1);
		}

		// starting base port
		{
			GRBLinExpr out_expr = 0;
			int port_i = param->portList[param->baseStartIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int out_index = param->getLegIndex(port_i, port_j);
				if (out_index != -1)
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
			}
			model.addConstr(out_expr == 1);
		}
		{
			// ending base port
			GRBLinExpr in_expr = 0;
			int port_i = param->portList[param->baseEndIndex].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				if (in_index != -1)
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
			}
			model.addConstr(in_expr == 1);
		}

		// optional ports
		for (int i = 0; i < param->addPortIndices.size(); i++)
		{
			GRBLinExpr in_expr = 0;
			GRBLinExpr out_expr = 0;
			int port_i = param->portList[param->addPortIndices[i]].portID;
			for (int j = 0; j < param->portNum; j++)
			{
				int port_j = param->portList[j].portID;
				if (port_i == port_j)
					continue;
				int in_index = param->getLegIndex(port_j, port_i);
				int out_index = param->getLegIndex(port_i, port_j);
				if (in_index != -1)
				{
					for (GRBVar& var : x_ij_s[in_index])
						in_expr += var;
				}
				if (out_index != -1)
				{
					for (GRBVar& var : x_ij_s[out_index])
						out_expr += var;
				}
			}
			model.addConstr(out_expr - y_i[i] == 0);
			model.addConstr(in_expr - y_i[i] == 0);
		}
	}

	// constraint (10)
	{
		for (int i = 0; i < param->prioList.size(); i++)
		{
			int id_i = param->prioList[i].first;
			int id_j = param->prioList[i].second;
			int leg_ind = param->getLegIndex(id_i, id_j);

			if (leg_ind != -1)
			{
				for (auto& var : x_ij_s[leg_ind])
				{
					model.addConstr(var - b_ij[i] <= 0);
				}
			}
		}
	}

	// constraint (11)
	{
		GRBLinExpr expr = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				int tail_ind = param->getPortIndex(param->legList[a].tailID);
				double voyTime = param->legList[a].confList[s].legTime + param->portList[tail_ind].portTime;
				expr += voyTime * x_ij_s[a][s];
			}
		}

		model.addConstr(expr <= param->maxVoyageTime);
	}

	// constraint (21)
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;
		for (int j = 0; j < param->portNum; j++)
		{
			if (j == param->baseStartIndex)
				continue;
			if (i == j)
				continue;

			int i_id = param->portList[i].portID;
			int j_id = param->portList[j].portID;

			GRBLinExpr expr1 = 0;
			int leg_ind1 = param->getLegIndex(i_id, j_id);
			if (leg_ind1 != -1)
			{
				for (int s = 0; s < param->legList[leg_ind1].confList.size(); s++)
				{
					double t_ij = param->legList[leg_ind1].confList[s].legTime + param->portList[i].portTime;
					expr1 += (param->maxVoyageTime + t_ij) * x_ij_s[leg_ind1][s];
				}
			}

			GRBLinExpr expr2 = 0;
			int leg_ind2 = param->getLegIndex(j_id, i_id);
			if (leg_ind2 != -1)
			{
				for (int s = 0; s < param->legList[leg_ind2].confList.size(); s++)
				{
					double t_ji = param->legList[leg_ind2].confList[s].legTime + param->portList[j].portTime;
					expr2 += (param->maxVoyageTime - t_ji) * x_ij_s[leg_ind2][s];
				}
			}

			model.addConstr(t_i[i] - t_i[j] + expr1 + expr2 <= param->maxVoyageTime);
		}
	}

	// constraint (22)
	GRBLinExpr expr1 = 0;
	for (int a = 0; a < param->legNum; a++)
	{
		int tail_id = param->legList[a].tailID;

		if (tail_id == param->baseStartID)
		{
			for (int s = 0; s < param->legList[a].confList.size(); s++)
			{
				double t_0j = param->legList[a].confList[s].legTime + param->portList[param->baseStartIndex].portTime;
				expr1 += t_0j * x_ij_s[a][s];
			}
		}
	}
	for (int i = 0; i < param->portNum; i++)
	{
		if (i == param->baseStartIndex)
			continue;

		int i_id = param->portList[i].portID;

		GRBLinExpr expr2 = 0;
		for (int a = 0; a < param->legNum; a++)
		{
			int tail_id = param->legList[a].tailID;
			int head_id = param->legList[a].headID;
			if (head_id == i_id && tail_id != param->baseStartID)
			{
				int j_ind = param->getPortIndex(tail_id);
				for (int s = 0; s < param->legList[a].confList.size(); s++)
				{
					double t_ji = param->legList[a].confList[s].legTime + param->portList[j_ind].portTime;
					expr2 += t_ji * x_ij_s[a][s];
				}
			}
		}
		model.addConstr(t_i[i] - expr1 - expr2 >= 0);
		model.addConstr(t_i[i] - param->maxVoyageTime <= 0);
	}

	// constraint (23)
	model.addConstr(t_i[param->baseStartIndex] == 0);

	// constraints (24)(25)
	for (int h = 0; h < param->scenarioList[sind].demandList.size(); h++)
	{
		double weight = param->scenarioList[sind].demandList[h].weight;
		double timeLimit = param->scenarioList[sind].demandList[h].maxTransTime;

		int o_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].originID);
		int d_ind = param->getPortIndex(param->scenarioList[sind].demandList[h].destinationID);

		model.addConstr(param->maxVoyageTime * (gama_h[h] - 1) - t_i[d_ind] + t_i[o_ind] <= 0);
		model.addConstr(t_i[d_ind] - t_i[o_ind] - timeLimit - param->maxVoyageTime * (1 - gama_h[h]) <= 0);
	}

	// cuts
	{
		vector<GRBLinExpr> y(param->portNum);
		for (int i = 0; i < param->portNum; i++)
		{
			y[i] = 0;
			int ind = -1;
			for (int j = 0; j < param->addPortIndices.size(); j++)
			{
				if (i == param->addPortIndices[j])
				{
					ind = j;
					break;
				}
			}
			if (ind != -1)
				y[i] += y_i[ind];
			else
				y[i] += 1;
		}
		for (int i = 0; i < param->portNum; i++)
		{
			if (i == param->baseStartIndex || i == param->baseEndIndex)
				continue;
			for (int j = 0; j < param->portNum; j++)
			{
				if (i == j)
					break;
				if (j == param->baseStartIndex || j == param->baseEndIndex)
					continue;

				int prioInd = param->getPrioIndex(param->portList[i].portID, param->portList[j].portID);

				if (prioInd != -1)
				{
					model.addConstr(y[i] - y[j] - b_ij[prioInd] <= 0);
					model.addConstr(b_ij[prioInd] - 1 - y[i] + y[j] <= 0);
					if (param->portList[i].portID < param->portList[j].portID)
						model.addConstr(1 - y[i] - y[j] - b_ij[prioInd] <= 0);
				}

			}
		}
	}

	// objective
	// objective function
	GRBLinExpr obj = 0;
	GRBLinExpr legCost = 0;
	GRBLinExpr portCost = 0;
	for (int i = 0; i < param->legNum; i++)
	{
		int tailID = param->legList[i].tailID;

		for (int s = 0; s < param->legList[i].confList.size(); s++)
		{
			double c_ij_s = param->legList[i].confList[s].legCost;
			legCost += c_ij_s * x_ij_s[i][s];
			portCost += param->portList[param->getPortIndex(tailID)].portCost*x_ij_s[i][s];
		}
	}

	obj = legCost + portCost + B;

	model.setObjective(obj, GRB_MINIMIZE);

	//model.getEnv().set(GRB_IntParam_OutputFlag, 0);
}

struct vars
{
	GRBVar1d b_ij;
	GRBVar2d x_ij_s;
	GRBVar1d y_i;
	GRBVar1d gama_h;
};

struct vals
{
	vector<int> b_ij;
	vector<vector<int>> x_ij_s;
	vector<int> y_i;
	vector<int> gama_h;
};

GRBLinExpr solve_subp(vars benders_var, vals val, int sind, double &objval)
{
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	int demandNum = val.gama_h.size();

	GRBVar1d f_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "f", param->legNum);
	GRBVar1d fv_ij = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "fv", param->legNum);
	GRBVar1d g_h = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "g", demandNum);

	vector<GRBConstr> x_constr = vector<GRBConstr>();
	vector<GRBConstr> y_constr = vector<GRBConstr>();
	vector<GRBConstr> b_constr = vector<GRBConstr>();
	vector<GRBConstr> gama_constr = vector<GRBConstr>();
	vector<GRBConstr> other_constr = vector<GRBConstr>();

	// constraint (12) (13)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			GRBLinExpr expr = 0;
			for (auto& var : val.x_ij_s[i])
				expr += var;
			x_constr.push_back(model.addConstr(f_ij[i] - param->shipWeightCapacity * expr <= 0));
			x_constr.push_back(model.addConstr(fv_ij[i] - param->shipVolumeCapacity * expr <= 0));
		}
	}

	// constraint (18) (19)
	{
		for (int i = 0; i < param->portNum; i++)
		{
			GRBLinExpr wlhs = 0;
			GRBLinExpr vlhs = 0;
			int port_i = param->portList[i].portID;
			for (int a = 0; a < param->legNum; a++)
			{
				int tailID = param->legList[a].tailID;
				int headID = param->legList[a].headID;

				if (headID == port_i)
				{
					wlhs += f_ij[a];
					vlhs += fv_ij[a];
				}

				if (tailID == port_i)
				{
					wlhs -= f_ij[a];
					vlhs -= fv_ij[a];
				}
			}

			GRBLinExpr wrhs = 0;
			GRBLinExpr vrhs = 0;
			for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
			{
				int originID = param->scenarioList[sind].demandList[k].originID;
				int destID = param->scenarioList[sind].demandList[k].destinationID;

				if (destID == port_i)
				{
					wrhs += g_h[k];
					vrhs += param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}

				if (originID == port_i)
				{
					wrhs -= g_h[k];
					vrhs -= param->scenarioList[sind].demandList[k].volumeRatio * g_h[k];
				}
			}
			model.addConstr(wlhs - wrhs == 0);
			model.addConstr(vlhs - vrhs == 0);
		}
	}

	// constraint (20)
	{
		for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
		{
			int originID = param->scenarioList[sind].demandList[k].originID;
			int destID = param->scenarioList[sind].demandList[k].destinationID;

			double demandAmount = param->scenarioList[sind].demandList[k].weight;

			if (originID == param->baseStartID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(destID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					y_constr.push_back(model.addConstr(g_h[k] - demandAmount * val.y_i[y_ind] <= 0));
				}
				else
				{
					other_constr.push_back(model.addConstr(g_h[k] - demandAmount <= 0));
				}
			}
			else if (destID == param->baseEndID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(originID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					y_constr.push_back(model.addConstr(g_h[k] - demandAmount * val.y_i[y_ind] <= 0));
				}
				else
				{
					other_constr.push_back(model.addConstr(g_h[k] - demandAmount <= 0));
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				b_constr.push_back(model.addConstr(g_h[k] - demandAmount * val.b_ij[prioInd] <= 0));

				int originInd = param->getPortIndex(originID);
				int destInd = param->getPortIndex(destID);

				int y_origin_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == originInd)
					{
						y_origin_ind = i;
						break;
					}
				}
				if (y_origin_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == originInd)
						{
							y_origin_ind = i;
							break;
						}
					}
					y_constr.push_back(model.addConstr(g_h[k] - demandAmount * val.y_i[y_origin_ind] <= 0));
				}

				int y_dest_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == destInd)
					{
						y_dest_ind = i;
						break;
					}
				}
				if (y_dest_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == destInd)
						{
							y_dest_ind = i;
							break;
						}
					}
					y_constr.push_back(model.addConstr(g_h[k] - demandAmount * val.y_i[y_dest_ind] <= 0));
				}
			}
		}
	}

	// constraints (24)(25)
	for (int h = 0; h < param->scenarioList[sind].demandList.size(); h++)
	{
		double weight = param->scenarioList[sind].demandList[h].weight;

		gama_constr.push_back(model.addConstr(g_h[h] - weight * val.gama_h[h] <= 0));
	}

	// objective
	GRBLinExpr obj = 0;
	GRBLinExpr frtRev = 0;
	for (int k = 0; k < demandNum; k++)
	{
		frtRev += param->scenarioList[sind].demandList[k].revenue*g_h[k];
	}
	obj = -frtRev;

	model.setObjective(obj, GRB_MINIMIZE);

	//model.write("sub.lp");

	model.getEnv().set(GRB_IntParam_OutputFlag, 0);

	model.optimize();

	vector<double> val_x_constr = getOptimalDualValue(x_constr);
	vector<double> val_y_constr = getOptimalDualValue(y_constr);
	vector<double> val_b_constr = getOptimalDualValue(b_constr);
	vector<double> val_gama_constr = getOptimalDualValue(gama_constr);
	vector<double> val_other_constr = getOptimalDualValue(other_constr);

	// generate expr
	GRBLinExpr benders_expr = 0;
	int x_iter = 0, y_iter = 0, gama_iter = 0, b_iter = 0, other_iter = 0;

	// constraint (12) (13)
	{
		for (int i = 0; i < param->legNum; i++)
		{
			GRBLinExpr expr = 0;
			for (auto& var : benders_var.x_ij_s[i])
				expr += var;
			//cout << std::to_string(expr.size()) << endl;
			benders_expr += val_x_constr[x_iter] * param->shipWeightCapacity*expr;
			x_iter++;
			benders_expr += val_x_constr[x_iter] * param->shipVolumeCapacity*expr;
			x_iter++;
		}
	}
	
	// constraint (20)
	{
		for (int k = 0; k < param->scenarioList[sind].demandList.size(); k++)
		{
			int originID = param->scenarioList[sind].demandList[k].originID;
			int destID = param->scenarioList[sind].demandList[k].destinationID;

			double demandAmount = param->scenarioList[sind].demandList[k].weight;

			if (originID == param->baseStartID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(destID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					benders_expr += val_y_constr[y_iter] * demandAmount*benders_var.y_i[y_ind];
					y_iter++;
				}
				else
				{
					benders_expr += demandAmount * val_other_constr[other_iter];
					other_iter++;
				}
			}
			else if (destID == param->baseEndID)
			{
				int y_ind = -1;
				int portInd = param->getPortIndex(originID);
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == portInd)
					{
						y_ind = i;
						break;
					}
				}
				if (y_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == portInd)
						{
							y_ind = i;
							break;
						}
					}
					benders_expr += val_y_constr[y_iter] * demandAmount*benders_var.y_i[y_ind];
					y_iter++;
				}
				else
				{
					benders_expr += demandAmount * val_other_constr[other_iter];
					other_iter++;
				}
			}
			else
			{
				int prioInd = param->getPrioIndex(originID, destID);
				benders_expr += val_b_constr[b_iter] * demandAmount*benders_var.b_ij[prioInd];
				b_iter++;

				int originInd = param->getPortIndex(originID);
				int destInd = param->getPortIndex(destID);

				int y_origin_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == originInd)
					{
						y_origin_ind = i;
						break;
					}
				}
				if (y_origin_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == originInd)
						{
							y_origin_ind = i;
							break;
						}
					}
					benders_expr += val_y_constr[y_iter] * demandAmount*benders_var.y_i[y_origin_ind];
					y_iter++;
					
				}

				int y_dest_ind = -1;
				for (int i = 0; i < param->fixPortIndices.size(); i++)
				{
					if (param->fixPortIndices[i] == destInd)
					{
						y_dest_ind = i;
						break;
					}
				}
				if (y_dest_ind == -1)
				{
					for (int i = 0; i < param->addPortIndices.size(); i++)
					{
						if (param->addPortIndices[i] == destInd)
						{
							y_dest_ind = i;
							break;
						}
					}
					benders_expr += val_y_constr[y_iter] * demandAmount*benders_var.y_i[y_dest_ind];
					y_iter++;
				}
			}
		}
	}

	for (int h = 0; h < param->scenarioList[sind].demandList.size(); h++)
	{
		double weight = param->scenarioList[sind].demandList[h].weight;

		//cout << h << endl;

		benders_expr += val_gama_constr[gama_iter] * weight*benders_var.gama_h[h];
		gama_iter++;

	}

	objval = model.get(GRB_DoubleAttr_ObjVal);

	return benders_expr;
}

MasterSolution SubLBProblem_benders::getSolution()
{
	double LB = -MAX_NUM;
	double UB = MAX_NUM;

	while(UB>LB+0.001)
	{
		model.optimize();

		LB = model.get(GRB_DoubleAttr_ObjVal);

		cout << "LB:" << LB << endl;

		double val_B = B.get(GRB_DoubleAttr_X);
		vals val;
		val.b_ij = getOptimalValue(b_ij);
		val.gama_h = getOptimalValue(gama_h);
		val.y_i = getOptimalValue(y_i);
		val.x_ij_s = getOptimalValue(x_ij_s);

		vars var;
		var.x_ij_s = x_ij_s;
		var.b_ij = b_ij;
		var.gama_h = gama_h;
		var.y_i = y_i;

		double obj_val = 0;
		GRBLinExpr expr = solve_subp(var, val, sind, obj_val);

		if(obj_val>val_B)
		{
			double t_UB = model.get(GRB_DoubleAttr_ObjVal) - val_B + obj_val;
			if (t_UB < UB)
				UB = t_UB;

			model.addConstr(B - expr >= 0);
		}

		cout << "UB:" << UB << endl;
	}

	MasterSolution ms;

	ms.a_ij = vector<int>(param->fixPrioList.size());
	for (int i = 0; i < ms.a_ij.size(); i++)
	{
		int id_i = param->fixPrioList[i].first;
		int id_j = param->fixPrioList[i].second;

		int prio2 = param->getPrioIndex(id_i, id_j);

		double val_a = b_ij[prio2].get(GRB_DoubleAttr_X);
		ms.a_ij[i] = Round(val_a);
	}
	ms.objVal = model.get(GRB_DoubleAttr_ObjVal);
	ms.objLB = ms.objVal;
	return ms;

}

void SubLBProblem_benders::setDualPenalty(vector<double>& omega)
{
	int demandNum = param->scenarioList[sind].demandList.size();

	// objective
	GRBLinExpr obj = 0;
	GRBLinExpr legCost = 0;
	GRBLinExpr portCost = 0;
	for (int i = 0; i < param->legNum; i++)
	{
		int headID = param->legList[i].headID;

		for (int s = 0; s < param->legList[i].confList.size(); s++)
		{
			double c_ij_s = param->legList[i].confList[s].legCost;
			legCost += c_ij_s * x_ij_s[i][s];
			portCost += param->portList[param->getPortIndex(headID)].portCost*x_ij_s[i][s];
		}
	}
	
	obj = legCost + portCost + B;

	for (int i = 0; i < param->fixPrioList.size(); i++)
	{
		int id_i = param->fixPrioList[i].first;
		int id_j = param->fixPrioList[i].second;

		int prio2 = param->getPrioIndex(id_i, id_j);
		obj += omega[i] * b_ij[prio2];
	}

	model.setObjective(obj, GRB_MINIMIZE);
}

void SubLBProblem_benders::removeSolution(MasterSolution& ms)
{
	GRBLinExpr expr = 0;
	for (int i = 0; i < param->fixPrioList.size(); i++)
	{
		int id_i = param->fixPrioList[i].first;
		int id_j = param->fixPrioList[i].second;

		int prio2 = param->getPrioIndex(id_i, id_j);

		if (ms.a_ij[i] == 1)
			expr += 1 - b_ij[prio2];
		else
			expr += b_ij[prio2];
	}
	model.addConstr(expr >= 2);

	feasible->removeSolution(ms);
}

void SubLBProblem_benders::removeSolutions(vector<MasterSolution>& vector_ms)
{
	for (int j = 0; j < vector_ms.size(); j++)
	{
		removeSolution(vector_ms[j]);
	}
}

